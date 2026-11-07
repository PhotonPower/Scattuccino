// Voxel transport, see cpp/include/scattuccino/core.hpp. Line-by-line mirror of voxels.py.
#include "scattuccino/core.hpp"

#include "scattuccino/photon.hpp"
#include "scattuccino/photon_pol.hpp"

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
  photon::PhaseView<double> v;
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
// plain views of the input for photon::trace (host memory), in the precision of the run
template <class Real>
struct Views {
  std::vector<photon::MaterialView<Real>> materials;
  std::vector<std::vector<Real>> store;  // converted tables of the tabulated phase functions
};

template <class Real>
Views<Real> make_views(const VoxelInput& in) {
  Views<Real> v;
  v.store.reserve(3 * in.materials.size());
  auto table = [&v](const std::vector<double>& src) -> const Real* {
    v.store.emplace_back(src.begin(), src.end());
    return v.store.back().data();
  };
  for (const auto& m : in.materials) {
    photon::MaterialView<Real> mv{};
    mv.mua = static_cast<Real>(m.mua);
    mv.mus = static_cast<Real>(m.mus);
    mv.n = static_cast<Real>(m.n);
    const double mt = m.mua + m.mus;
    mv.absfrac = mt > 0 ? m.mua / mt : 0.0;
    mv.phase.kind = static_cast<int>(m.phase.kind);
    mv.phase.p1 = static_cast<Real>(m.phase.p1);
    mv.phase.p2 = static_cast<Real>(m.phase.p2);
    mv.phase.mu = table(m.phase.mu);
    mv.phase.p = table(m.phase.p);
    mv.phase.cdf = table(m.phase.cdf);
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

namespace {
template <class Real>
VoxelOutput run(const VoxelInput& in) {
  const Setup st = prepare(in);
  const Views<Real> views = make_views<Real>(in);
  const std::size_t nm = in.materials.size();
  const std::size_t nvox = static_cast<std::size_t>(in.nx) * in.ny * in.nz;
  photon::Geometry<Real> geo{in.nx, in.ny, in.nz,
                             {static_cast<Real>(in.dx), static_cast<Real>(in.dy), static_cast<Real>(in.dz)},
                             {static_cast<Real>(st.lo[0]), static_cast<Real>(st.lo[1]), 0},
                             in.labels.data(), views.materials.data(),
                             static_cast<Real>(in.n_above), static_cast<Real>(in.n_below),
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
      const photon::Exit e = photon::trace<Real>(geo, in.seed, static_cast<std::uint64_t>(id), st.sin_t0,
                                           st.cos_t0, 1.0 - st.r_specular, st.i0, st.j0, deposit, add_path);
      if (e.side == 2) {  // roulette
        ++o.n_roulette;
        continue;
      }
      ++o.n_exit;
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
    out.n_exit += o.n_exit;
    out.n_roulette += o.n_roulette;
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

}  // namespace

VoxelOutput simulate_voxels(const VoxelInput& in) {
  switch (in.precision) {
    case 0:
      return run<double>(in);
    case 1:
      return run<float>(in);
    default:
      throw std::invalid_argument("precision must be 0 (double) or 1 (single)");
  }
}

PolOutput simulate_voxels_polarised(const PolInput& pin) {
  const VoxelInput& in = pin.voxels;
  const Setup st = prepare(in);  // checks labels, entry voxel
  const std::size_t nm = in.materials.size();
  if (pin.mueller.size() != nm) throw std::invalid_argument("one Mueller table per material");
  std::vector<photon::PolMaterialView> mats;
  for (std::size_t m = 0; m < nm; ++m) {
    const Material& M = in.materials[m];
    const MuellerTable& T = pin.mueller[m];
    if (M.phase.kind != PhaseKind::Tabulated) throw std::invalid_argument("polarised: phase must be tabulated S11");
    photon::PolMaterialView mv{M.mua, M.mus, M.n, {}, {}};
    mv.phase.kind = 3;
    mv.phase.mu = M.phase.mu.data();
    mv.phase.p = M.phase.p.data();
    mv.phase.cdf = M.phase.cdf.data();
    mv.phase.n = static_cast<int>(M.phase.mu.size());
    mv.mueller = {T.mu.data(), T.m12.data(), T.m33.data(), T.m34.data(), static_cast<int>(T.mu.size()), T.pure};
    mats.push_back(mv);
  }
  const photon::PolGeometry geo{in.nx, in.ny, in.nz, {in.dx, in.dy, in.dz}, {st.lo[0], st.lo[1], 0.0},
                                in.labels.data(), mats.data(), in.n_above, in.n_below, in.periodic,
                                in.w_threshold, in.roulette_m};
  // polarised entry, as simulate_polarised (polarisation.py)
  const double n_slab = in.materials[static_cast<std::size_t>(in.labels[(static_cast<std::size_t>(st.i0) * in.ny + st.j0) * in.nz])].n;
  const double ci = in.cos_incident, si = std::sqrt(1 - ci * ci);
  const double u0[3] = {si, 0.0, ci};
  double v0[3] = {1.0, 0.0, 0.0};
  if (si > 0) v0[0] = ci, v0[2] = -si;
  double s_in[4];
  for (int k = 0; k < 4; ++k) s_in[k] = pin.stokes_in[static_cast<std::size_t>(k)] / pin.stokes_in[0];
  const double zhat[3] = {0.0, 0.0, 1.0};
  double p_hat[3], s_hat[3], stp[4];
  photon::frame_to_plane(u0, v0, s_in, zhat, p_hat, s_hat, stp);
  photon::cplx rp, rs, tp, ts, ct;
  photon::fresnel_amplitudes(in.n_above, n_slab, ci, rp, rs, tp, ts, ct);
  double s_r[4], s_t[4];
  photon::apply_amplitudes(rp, rs, stp, s_r);
  photon::apply_amplitudes(tp, ts, stp, s_t);
  const double big_r = s_r[0];
  const double ratio = in.n_above / n_slab;
  double ut[3];
  for (int a = 0; a < 3; ++a) ut[a] = ratio * u0[a] + (ct.real() - ratio * ci) * zhat[a];
  const double nu = photon::norm3(ut);
  for (double& x : ut) x /= nu;
  double vt[3];
  photon::cross(s_hat, ut, vt);
  double S0[4];
  for (int k = 0; k < 4; ++k) S0[k] = s_t[k] / s_t[0];

  PolOutput out;
  out.specular_r = big_r;
  for (int k = 0; k < 4; ++k) out.specular_stokes.push_back(big_r > 0 ? s_r[k] / s_r[0] : s_r[k]);
  out.absorbed.assign(nm, 0.0);
  const int threads = std::max(in.threads, 1);
  std::vector<PolOutput> parts(static_cast<std::size_t>(threads));
#pragma omp parallel for num_threads(threads) schedule(static)
  for (int t = 0; t < threads; ++t) {
    PolOutput& o = parts[static_cast<std::size_t>(t)];
    o.absorbed.assign(nm, 0.0);
    const std::int64_t begin = in.packets * t / threads, end = in.packets * (t + 1) / threads;
    for (std::int64_t id = begin; id < end; ++id) {
      auto deposit = [&](int mat, int, int, int, double dw) { o.absorbed[static_cast<std::size_t>(mat)] += dw; };
      const photon::PolExit e = photon::trace_pol(geo, in.seed, static_cast<std::uint64_t>(id), ut, vt, S0,
                                                  1.0 - big_r, st.i0, st.j0, deposit);
      if (e.side == 2 || !in.record) continue;
      o.id.push_back(id);
      o.side.push_back(e.side);
      o.x.push_back(e.x);
      o.y.push_back(e.y);
      o.u.insert(o.u.end(), e.u, e.u + 3);
      o.v.insert(o.v.end(), e.v, e.v + 3);
      o.stokes.insert(o.stokes.end(), e.S, e.S + 4);
      o.w.push_back(e.w);
      o.scatterings.push_back(e.scatterings);
    }
  }
  for (const auto& o : parts) {
    for (std::size_t k = 0; k < nm; ++k) out.absorbed[k] += o.absorbed[k];
    out.id.insert(out.id.end(), o.id.begin(), o.id.end());
    out.side.insert(out.side.end(), o.side.begin(), o.side.end());
    for (auto [dst, src] : {std::pair{&out.x, &o.x}, {&out.y, &o.y}, {&out.u, &o.u}, {&out.v, &o.v},
                            {&out.stokes, &o.stokes}, {&out.w, &o.w}})
      dst->insert(dst->end(), src->begin(), src->end());
    out.scatterings.insert(out.scatterings.end(), o.scatterings.begin(), o.scatterings.end());
  }
  return out;
}

}  // namespace scattuccino
