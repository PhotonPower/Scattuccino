// CUDA backend (stage S5b): one thread per packet, running photon::trace (the same code as the
// CPU core) with the packet's own random stream; accumulators with atomicAdd. Exits are written
// to per-packet slots (each packet leaves at most once) and compacted on the host.
#include <cuda_runtime.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <memory>
#include <vector>

#include "scattuccino/core.hpp"
#include "scattuccino/photon.hpp"

namespace scattuccino {

namespace {

void check(cudaError_t e, const char* what) {
  if (e != cudaSuccess) throw std::runtime_error(std::string("CUDA error in ") + what + ": " + cudaGetErrorString(e));
}

template <class T>
struct Buffer {
  T* ptr = nullptr;
  std::size_t n = 0;
  explicit Buffer(std::size_t count) : n(count) {
    if (n) check(cudaMalloc(&ptr, n * sizeof(T)), "cudaMalloc");
  }
  ~Buffer() {
    if (ptr) cudaFree(ptr);
  }
  Buffer(const Buffer&) = delete;
  Buffer& operator=(const Buffer&) = delete;
  void upload(const T* h) { if (n) check(cudaMemcpy(ptr, h, n * sizeof(T), cudaMemcpyHostToDevice), "upload"); }
  void download(T* h) const { if (n) check(cudaMemcpy(h, ptr, n * sizeof(T), cudaMemcpyDeviceToHost), "download"); }
  void zero() { if (n) check(cudaMemset(ptr, 0, n * sizeof(T)), "memset"); }
};

struct DeviceExit {
  int side;
  double x, y, ux, uy, uz, w, zmax;
  long long scatterings;
};

constexpr int kLocalMaterials = 16;  // per-thread accumulators; more materials use atomics directly

// Persistent threads: every thread takes the next packet from a global counter until all are
// done (no idle warps at the end of a launch with long-lived packets). Absorption per material is
// summed in registers and flushed once per packet; the voxel grid is scored only on request.
template <class Real, class Pos>
__global__ void k_trace(photon::Geometry<Real> geo, std::uint64_t seed, long long first, long long count,
                        unsigned long long* next, double sin_t0, double cos_t0, double w0, int i0,
                        int j0, int nm, bool score_voxels, double* absorbed, double* absorbed_voxels,
                        double* sums, unsigned long long* counts, bool record, DeviceExit* exits,
                        double* paths) {
  double acc[kLocalMaterials];
  double s_r = 0.0, s_r2 = 0.0, s_t = 0.0, s_t2 = 0.0, s_side = 0.0;
  unsigned long long n_exit = 0, n_roulette = 0;
  while (true) {
    const long long t = static_cast<long long>(atomicAdd(next, 1ULL));
    if (t >= count) break;
    const long long id = first + t;
    for (int k = 0; k < kLocalMaterials; ++k) acc[k] = 0.0;
    double* my_path = record ? paths + t * nm : nullptr;
    auto deposit = [&](int mat, int i, int j, int k, double dw) {
      if (mat < kLocalMaterials) acc[mat] += dw; else atomicAdd(&absorbed[mat], dw);
      if (score_voxels) atomicAdd(&absorbed_voxels[(static_cast<long long>(i) * geo.ny + j) * geo.nz + k], dw);
    };
    auto add_path = [&](int mat, double d) {
      if (my_path) my_path[mat] += d;
    };
    const photon::Exit e = photon::trace<Real, Pos>(geo, seed, static_cast<std::uint64_t>(id), sin_t0, cos_t0, w0, i0,
                                         j0, deposit, add_path);
    for (int k = 0; k < kLocalMaterials && k < nm; ++k)
      if (acc[k] != 0.0) atomicAdd(&absorbed[k], acc[k]);
    if (e.side == 2) ++n_roulette; else ++n_exit;
    if (e.side == -1) { s_r += e.w; s_r2 += e.w * e.w; }
    if (e.side == 1) { s_t += e.w; s_t2 += e.w * e.w; }
    if (e.side == 0) s_side += e.w;
    if (record) exits[t] = DeviceExit{e.side, e.x, e.y, e.ux, e.uy, e.uz, e.w, e.zmax, e.scatterings};
  }
  atomicAdd(&sums[0], s_r);
  atomicAdd(&sums[1], s_r2);
  atomicAdd(&sums[2], s_t);
  atomicAdd(&sums[3], s_t2);
  atomicAdd(&sums[4], s_side);
  atomicAdd(&counts[0], n_exit);
  atomicAdd(&counts[1], n_roulette);
}

}  // namespace

int cuda_device_count() {
  int n = 0;
  if (cudaGetDeviceCount(&n) != cudaSuccess) {
    cudaGetLastError();
    return 0;
  }
  return n;
}

namespace {

template <class Real, class Pos>
VoxelOutput run_cuda(const VoxelInput& in, int device) {
  check(cudaSetDevice(device), "cudaSetDevice");
  const Setup st = prepare(in);
  const std::size_t nm = in.materials.size();
  const std::size_t nvox = static_cast<std::size_t>(in.nx) * in.ny * in.nz;
  // device copies: labels, tabulated phase arrays, materials
  Buffer<std::int32_t> d_labels(nvox);
  d_labels.upload(in.labels.data());
  std::vector<photon::MaterialView<Real>> mats(nm);
  std::vector<std::unique_ptr<Buffer<Real>>> tables;
  for (std::size_t m = 0; m < nm; ++m) {
    const Material& M = in.materials[m];
    photon::MaterialView<Real> mv{};
    mv.mua = static_cast<Real>(M.mua);
    mv.mus = static_cast<Real>(M.mus);
    mv.n = static_cast<Real>(M.n);
    const double mt = M.mua + M.mus;
    mv.absfrac = mt > 0 ? M.mua / mt : 0.0;
    mv.phase.kind = static_cast<int>(M.phase.kind);
    mv.phase.p1 = static_cast<Real>(M.phase.p1);
    mv.phase.p2 = static_cast<Real>(M.phase.p2);
    if (M.phase.kind == PhaseKind::Tabulated) {
      const std::vector<double>* arrays[3] = {&M.phase.mu, &M.phase.p, &M.phase.cdf};
      const Real* ptrs[3];
      for (int a = 0; a < 3; ++a) {
        const std::vector<Real> conv(arrays[a]->begin(), arrays[a]->end());
        tables.push_back(std::make_unique<Buffer<Real>>(conv.size()));
        tables.back()->upload(conv.data());
        ptrs[a] = tables.back()->ptr;
      }
      mv.phase.mu = ptrs[0];
      mv.phase.p = ptrs[1];
      mv.phase.cdf = ptrs[2];
      mv.phase.n = static_cast<int>(M.phase.mu.size());
    }
    mats[m] = mv;
  }
  Buffer<photon::MaterialView<Real>> d_mats(nm);
  d_mats.upload(mats.data());
  photon::Geometry<Real> geo{in.nx, in.ny, in.nz,
                             {static_cast<Real>(in.dx), static_cast<Real>(in.dy), static_cast<Real>(in.dz)},
                             {static_cast<Real>(st.lo[0]), static_cast<Real>(st.lo[1]), 0},
                             d_labels.ptr, d_mats.ptr,
                             static_cast<Real>(in.n_above), static_cast<Real>(in.n_below),
                             in.periodic, in.w_threshold, in.roulette_m};
  Buffer<double> d_abs(nm), d_vox(nvox), d_sums(5);
  Buffer<unsigned long long> d_next(1), d_counts(2);
  int sms = 0;
  check(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device), "attribute");
  d_abs.zero();
  d_vox.zero();
  d_sums.zero();
  d_counts.zero();
  VoxelOutput out;
  out.r_specular = st.r_specular;
  const long long chunk = in.record ? 1 << 20 : 1 << 24;  // packets per launch
  std::unique_ptr<Buffer<DeviceExit>> d_exits;
  std::unique_ptr<Buffer<double>> d_paths;
  if (in.record) {
    d_exits = std::make_unique<Buffer<DeviceExit>>(static_cast<std::size_t>(chunk));
    d_paths = std::make_unique<Buffer<double>>(static_cast<std::size_t>(chunk) * nm);
  }
  std::vector<DeviceExit> h_exits;
  std::vector<double> h_paths;
  for (long long first = 0; first < in.packets; first += chunk) {
    const long long count = std::min<long long>(chunk, in.packets - first);
    if (in.record) d_paths->zero();
    d_next.zero();
    const int block = 128;
    const long long wanted = (count + block - 1) / block;
    const unsigned grid = static_cast<unsigned>(std::min<long long>(wanted, static_cast<long long>(sms) * 16));
    k_trace<Real, Pos><<<grid, block>>>(geo, in.seed, first, count, d_next.ptr, st.sin_t0, st.cos_t0, 1.0 - st.r_specular,
                             st.i0, st.j0, static_cast<int>(nm), in.score_voxels, d_abs.ptr, d_vox.ptr,
                             d_sums.ptr, d_counts.ptr, in.record, in.record ? d_exits->ptr : nullptr,
                             in.record ? d_paths->ptr : nullptr);
    check(cudaGetLastError(), "k_trace");
    check(cudaDeviceSynchronize(), "synchronize");
    if (in.record) {
      h_exits.resize(static_cast<std::size_t>(count));
      h_paths.resize(static_cast<std::size_t>(count) * nm);
      check(cudaMemcpy(h_exits.data(), d_exits->ptr, h_exits.size() * sizeof(DeviceExit), cudaMemcpyDeviceToHost),
            "download exits");
      check(cudaMemcpy(h_paths.data(), d_paths->ptr, h_paths.size() * sizeof(double), cudaMemcpyDeviceToHost),
            "download paths");
      for (long long t = 0; t < count; ++t) {
        const DeviceExit& e = h_exits[static_cast<std::size_t>(t)];
        if (e.side == 2) continue;
        out.id.push_back(first + t);
        out.side.push_back(e.side);
        out.x.push_back(e.x);
        out.y.push_back(e.y);
        out.ux.push_back(e.ux);
        out.uy.push_back(e.uy);
        out.uz.push_back(e.uz);
        out.w.push_back(e.w);
        out.zmax.push_back(e.zmax);
        out.path.insert(out.path.end(), h_paths.begin() + t * static_cast<long long>(nm),
                        h_paths.begin() + (t + 1) * static_cast<long long>(nm));
        out.scatterings.push_back(e.scatterings);
      }
    }
  }
  out.absorbed.resize(nm);
  out.absorbed_voxels.resize(nvox);
  d_abs.download(out.absorbed.data());
  d_vox.download(out.absorbed_voxels.data());
  double sums[5];
  d_sums.download(sums);
  out.sum_r = sums[0];
  out.sum_r2 = sums[1];
  out.sum_t = sums[2];
  out.sum_t2 = sums[3];
  out.sum_side = sums[4];
  unsigned long long counts[2];
  d_counts.download(counts);
  out.n_exit = static_cast<std::int64_t>(counts[0]);
  out.n_roulette = static_cast<std::int64_t>(counts[1]);
  return out;
}

}  // namespace

VoxelOutput simulate_voxels_cuda(const VoxelInput& in, int device) {
  switch (in.precision) {
    case 0:
      return run_cuda<double, double>(in, device);
    case 1:
      return run_cuda<float, float>(in, device);
    case 2:
      return run_cuda<float, double>(in, device);
    default:
      throw std::invalid_argument("precision must be 0 (double), 1 or 2 (single)");
  }
}

}  // namespace scattuccino
