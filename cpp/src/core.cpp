// Voxel transport, see cpp/include/scattuccino/core.hpp. Line-by-line mirror of voxels.py.
#include "scattuccino/core.hpp"

#include "scattuccino/photon.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace scattuccino {


CounterRandom::CounterRandom(std::uint64_t seed, std::uint64_t id)
    : key(photon::Rng(seed, id).key) {}

double CounterRandom::uniform() {
  photon::Rng r(0, 0);
  r.key = key;
  r.counter = counter++;
  return r.uniform();
}

double Phase::sample(double xi) const {
  photon::PhaseView v;
  v.kind = static_cast<int>(kind);
  v.p1 = p1;
  v.p2 = p2;
  v.mu = mu.data();
  v.p = p.data();
  v.cdf = cdf.data();
  v.n = static_cast<int>(mu.size());
  return photon::sample(v, xi);
}

namespace {
// plain views of the input for photon::trace (host memory)
struct Views {
  std::vector<photon::MaterialView> materials;
};

Views make_views(const VoxelInput& in) {
  Views v;
  for (const auto& m : in.materials) {
    photon::MaterialView mv{m.mua, m.mus, m.n, {}};
    mv.phase.kind = static_cast<int>(m.phase.kind);
    mv.phase.p1 = m.phase.p1;
    mv.phase.p2 = m.phase.p2;
    mv.phase.mu = m.phase.mu.data();
    mv.phase.p = m.phase.p.data();
    mv.phase.cdf = m.phase.cdf.data();
    mv.phase.n = static_cast<int>(m.phase.mu.size());
    v.materials.push_back(mv);
  }
  return v;
}
}  // namespace

Setup prepare(const VoxelInput& in) {
  const std::size_t nvox = static_cast<std::size_t>(in.nx) * in.ny * in.nz;
  if (in.labels.size() != nvox) throw std::invalid_argument("labels must have nx·ny·nz entries");
  for (auto l : in.labels)
    if (l < 0 || static_cast<std::size_t>(l) >= in.materials.size())
      throw std::invalid_argument("label out of range");
  for (const auto& m : in.materials)
    if (m.phase.kind == PhaseKind::Tabulated && (m.phase.mu.size() < 2 || m.phase.mu.size() != m.phase.p.size() ||
                                                 m.phase.mu.size() != m.phase.cdf.size()))
      throw std::invalid_argument("tabulated phase function needs mu, p, cdf of equal length ≥ 2");
  Setup s;
  s.lo[0] = -in.nx * in.dx / 2;
  s.lo[1] = -in.ny * in.dy / 2;
  s.i0 = std::min(static_cast<int>(-s.lo[0] / in.dx), in.nx - 1);
  s.j0 = std::min(static_cast<int>(-s.lo[1] / in.dy), in.ny - 1);
  const Material& entry =
      in.materials[static_cast<std::size_t>(in.labels[(static_cast<std::size_t>(s.i0) * in.ny + s.j0) * in.nz])];
  double r_sp, cos_t0;
  photon::fresnel(in.n_above, entry.n, in.cos_incident, r_sp, cos_t0);
  s.r_specular = r_sp;
  s.cos_t0 = cos_t0;
  s.sin_t0 = in.n_above / entry.n * std::sqrt(1 - in.cos_incident * in.cos_incident);
  return s;
}

VoxelOutput simulate_voxels(const VoxelInput& in) {
  const Setup st = prepare(in);
  const Views views = make_views(in);
  const std::size_t nm = in.materials.size();
  const std::size_t nvox = static_cast<std::size_t>(in.nx) * in.ny * in.nz;
  photon::Geometry geo{in.nx, in.ny, in.nz, {in.dx, in.dy, in.dz}, {st.lo[0], st.lo[1], 0.0},
                       in.labels.data(), views.materials.data(), in.n_above, in.n_below,
                       in.periodic, in.w_threshold, in.roulette_m};
  VoxelOutput out;
  out.r_specular = st.r_specular;
  out.absorbed.assign(nm, 0.0);
  out.absorbed_voxels.assign(nvox, 0.0);
  const int threads = std::max(in.threads, 1);
  std::vector<VoxelOutput> parts(static_cast<std::size_t>(threads));
#pragma omp parallel for num_threads(threads) schedule(static)
  for (int t = 0; t < threads; ++t) {
    VoxelOutput& o = parts[static_cast<std::size_t>(t)];
    o.absorbed.assign(nm, 0.0);
    o.absorbed_voxels.assign(nvox, 0.0);
    const std::int64_t begin = in.packets * t / threads, end = in.packets * (t + 1) / threads;
    std::vector<double> path(nm);
    for (std::int64_t id = begin; id < end; ++id) {
      std::fill(path.begin(), path.end(), 0.0);
      auto deposit = [&](int mat, int i, int j, int k, double dw) {
        o.absorbed[static_cast<std::size_t>(mat)] += dw;
        o.absorbed_voxels[(static_cast<std::size_t>(i) * in.ny + j) * in.nz + k] += dw;
      };
      auto add_path = [&](int mat, double d) { path[static_cast<std::size_t>(mat)] += d; };
      const photon::Exit e = photon::trace(geo, in.seed, static_cast<std::uint64_t>(id), st.sin_t0,
                                           st.cos_t0, 1.0 - st.r_specular, st.i0, st.j0, deposit, add_path);
      if (e.side == 2) continue;  // roulette
      const double w2 = e.w * e.w;
      if (e.side == -1) o.sum_r += e.w, o.sum_r2 += w2;
      if (e.side == 1) o.sum_t += e.w, o.sum_t2 += w2;
      if (e.side == 0) o.sum_side += e.w;
      if (in.record) {
        o.id.push_back(id);
        o.side.push_back(e.side);
        o.x.push_back(e.x);
        o.y.push_back(e.y);
        o.ux.push_back(e.ux);
        o.uy.push_back(e.uy);
        o.uz.push_back(e.uz);
        o.w.push_back(e.w);
        o.zmax.push_back(e.zmax);
        o.path.insert(o.path.end(), path.begin(), path.end());
        o.scatterings.push_back(e.scatterings);
      }
    }
  }
  for (const auto& o : parts) {
    for (std::size_t k = 0; k < nm; ++k) out.absorbed[k] += o.absorbed[k];
    for (std::size_t k = 0; k < nvox; ++k) out.absorbed_voxels[k] += o.absorbed_voxels[k];
    out.sum_r += o.sum_r;
    out.sum_r2 += o.sum_r2;
    out.sum_t += o.sum_t;
    out.sum_t2 += o.sum_t2;
    out.sum_side += o.sum_side;
    out.id.insert(out.id.end(), o.id.begin(), o.id.end());
    out.side.insert(out.side.end(), o.side.begin(), o.side.end());
    out.x.insert(out.x.end(), o.x.begin(), o.x.end());
    out.y.insert(out.y.end(), o.y.begin(), o.y.end());
    out.ux.insert(out.ux.end(), o.ux.begin(), o.ux.end());
    out.uy.insert(out.uy.end(), o.uy.begin(), o.uy.end());
    out.uz.insert(out.uz.end(), o.uz.begin(), o.uz.end());
    out.w.insert(out.w.end(), o.w.begin(), o.w.end());
    out.zmax.insert(out.zmax.end(), o.zmax.begin(), o.zmax.end());
    out.path.insert(out.path.end(), o.path.begin(), o.path.end());
    out.scatterings.insert(out.scatterings.end(), o.scatterings.begin(), o.scatterings.end());
  }
  return out;
}

}  // namespace scattuccino
