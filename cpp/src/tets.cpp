// Tetrahedral transport, see cpp/include/scattuccino/tets.hpp and photon_tet.hpp.
#include "scattuccino/tets.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>

#include "scattuccino/photon_tet.hpp"

namespace scattuccino {

namespace {
using V3 = std::array<double, 3>;
V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
V3 cross(const V3& a, const V3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
}  // namespace

TetOutput simulate_tets(const TetInput& in) {
  const std::size_t nt = in.tets.size() / 4, nn = in.nodes.size() / 3;
  if (in.tets.size() != nt * 4 || in.nodes.size() != nn * 3 || in.labels.size() != nt)
    throw std::invalid_argument("nodes (N, 3), tets (M, 4), labels (M) expected");
  const std::size_t nm = in.materials.size();
  for (auto l : in.labels)
    if (l < 0 || static_cast<std::size_t>(l) >= nm) throw std::invalid_argument("label out of range");
  auto node = [&](int i) { return V3{in.nodes[3 * static_cast<std::size_t>(i)], in.nodes[3 * static_cast<std::size_t>(i) + 1],
                                     in.nodes[3 * static_cast<std::size_t>(i) + 2]}; };
  // faces opposite vertex f: the other three vertices
  constexpr int opp[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  std::vector<double> normal(nt * 12), offset(nt * 4), ext(nt * 4, 1.0);
  std::vector<int> neighbour(nt * 4, -1);
  std::map<std::array<int, 3>, std::pair<int, int>> faces;
  for (std::size_t t = 0; t < nt; ++t) {
    const int* v = &in.tets[4 * t];
    for (int f = 0; f < 4; ++f) {
      const V3 a = node(v[opp[f][0]]), b = node(v[opp[f][1]]), c = node(v[opp[f][2]]), d = node(v[f]);
      V3 n = cross(sub(b, a), sub(c, a));
      const double len = std::sqrt(dot(n, n));
      if (len == 0) throw std::invalid_argument("degenerate tetrahedron");
      for (auto& x : n) x /= len;
      if (dot(n, sub(d, a)) > 0)
        for (auto& x : n) x = -x;  // outward: away from the opposite vertex
      for (int k = 0; k < 3; ++k) normal[t * 12 + static_cast<std::size_t>(f) * 3 + static_cast<std::size_t>(k)] = n[static_cast<std::size_t>(k)];
      offset[t * 4 + static_cast<std::size_t>(f)] = dot(n, a);
      std::array<int, 3> key{v[opp[f][0]], v[opp[f][1]], v[opp[f][2]]};
      std::sort(key.begin(), key.end());
      auto it = faces.find(key);
      if (it == faces.end()) {
        faces.emplace(key, std::make_pair(static_cast<int>(t), f));
      } else {
        const auto [t2, f2] = it->second;
        neighbour[t * 4 + static_cast<std::size_t>(f)] = t2;
        neighbour[static_cast<std::size_t>(t2) * 4 + static_cast<std::size_t>(f2)] = static_cast<int>(t);
        faces.erase(it);
      }
    }
  }
  for (std::size_t q = 0; q < nt * 4; ++q) {
    if (neighbour[q] >= 0) continue;
    const double nz = normal[q * 3 + 2];
    ext[q] = in.n_all > 0 ? in.n_all : (nz < -0.5 ? in.n_top : (nz > 0.5 ? in.n_bottom : in.n_side));
  }
  // entry: boundary face containing the entry point, entered against its outward normal
  const V3 p{in.entry_point[0], in.entry_point[1], in.entry_point[2]};
  V3 d{in.direction[0], in.direction[1], in.direction[2]};
  {
    const double l = std::sqrt(dot(d, d));
    for (auto& x : d) x /= l;
  }
  int tet0 = -1, face0 = -1;
  double best = HUGE_VAL;
  for (std::size_t t = 0; t < nt && tet0 < 0; ++t)
    for (int f = 0; f < 4; ++f) {
      const std::size_t q = t * 4 + static_cast<std::size_t>(f);
      if (neighbour[q] >= 0) continue;
      const V3 n{normal[q * 3], normal[q * 3 + 1], normal[q * 3 + 2]};
      if (dot(n, d) >= 0) continue;
      const double dist = std::abs(dot(n, p) - offset[q]);
      if (dist > 1e-9) continue;
      // inside the triangle: barycentric test
      const int* v = &in.tets[4 * t];
      const V3 a = node(v[opp[f][0]]), b = node(v[opp[f][1]]), c = node(v[opp[f][2]]);
      const V3 c1 = cross(sub(b, a), sub(p, a)), c2 = cross(sub(c, b), sub(p, b)), c3 = cross(sub(a, c), sub(p, c));
      if (dot(c1, n) >= -1e-12 && dot(c2, n) >= -1e-12 && dot(c3, n) >= -1e-12 && dist < best) {
        best = dist;
        tet0 = static_cast<int>(t);
        face0 = f;
      }
    }
  if (tet0 < 0) throw std::invalid_argument("entry point is not on a boundary face facing the beam");
  const std::size_t q0 = static_cast<std::size_t>(tet0) * 4 + static_cast<std::size_t>(face0);
  double nin[3] = {-normal[q0 * 3], -normal[q0 * 3 + 1], -normal[q0 * 3 + 2]};  // into the mesh
  const double n_out = ext[q0], n_in = in.materials[static_cast<std::size_t>(in.labels[static_cast<std::size_t>(tet0)])].n;
  const double cos_i = d[0] * nin[0] + d[1] * nin[1] + d[2] * nin[2];
  double r_sp, cos_t;
  photon::fresnel(n_out, n_in, cos_i, r_sp, cos_t);
  double u0[3] = {d[0], d[1], d[2]};
  photon::refract(u0, nin, cos_i, cos_t, n_out / n_in);

  std::vector<photon::MaterialView> mats;
  for (const auto& M : in.materials) {
    photon::MaterialView mv{M.mua, M.mus, M.n, {}};
    mv.phase.kind = static_cast<int>(M.phase.kind);
    mv.phase.p1 = M.phase.p1;
    mv.phase.p2 = M.phase.p2;
    mv.phase.mu = M.phase.mu.data();
    mv.phase.p = M.phase.p.data();
    mv.phase.cdf = M.phase.cdf.data();
    mv.phase.n = static_cast<int>(M.phase.mu.size());
    mats.push_back(mv);
  }
  const photon::TetMesh mesh{normal.data(), offset.data(), neighbour.data(), ext.data(), in.labels.data(),
                             mats.data(), in.w_threshold, in.roulette_m};
  TetOutput out;
  out.r_specular = r_sp;
  out.entry_tet = tet0;
  out.absorbed.assign(nm, 0.0);
  out.absorbed_tets.assign(nt, 0.0);
  const int threads = std::max(in.threads, 1);
  std::vector<TetOutput> parts(static_cast<std::size_t>(threads));
#pragma omp parallel for num_threads(threads) schedule(static)
  for (int t = 0; t < threads; ++t) {
    TetOutput& o = parts[static_cast<std::size_t>(t)];
    o.absorbed.assign(nm, 0.0);
    o.absorbed_tets.assign(nt, 0.0);
    std::vector<double> path(nm);
    const std::int64_t begin = in.packets * t / threads, end = in.packets * (t + 1) / threads;
    for (std::int64_t id = begin; id < end; ++id) {
      std::fill(path.begin(), path.end(), 0.0);
      auto deposit = [&](int mat, int tet, double dw) {
        o.absorbed[static_cast<std::size_t>(mat)] += dw;
        o.absorbed_tets[static_cast<std::size_t>(tet)] += dw;
      };
      auto add_path = [&](int mat, double dd) { path[static_cast<std::size_t>(mat)] += dd; };
      const photon::TetExit e = photon::trace_tet(mesh, in.seed, static_cast<std::uint64_t>(id), in.entry_point.data(),
                                                  u0, 1.0 - r_sp, tet0, deposit, add_path);
      if (e.side == 2 || !in.record) continue;
      o.id.push_back(id);
      o.side.push_back(e.side);
      o.x.push_back(e.x);
      o.y.push_back(e.y);
      o.z.push_back(e.z);
      o.ux.push_back(e.ux);
      o.uy.push_back(e.uy);
      o.uz.push_back(e.uz);
      o.w.push_back(e.w);
      o.path.insert(o.path.end(), path.begin(), path.end());
      o.scatterings.push_back(e.scatterings);
    }
  }
  for (const auto& o : parts) {
    for (std::size_t k = 0; k < nm; ++k) out.absorbed[k] += o.absorbed[k];
    for (std::size_t k = 0; k < nt; ++k) out.absorbed_tets[k] += o.absorbed_tets[k];
    out.id.insert(out.id.end(), o.id.begin(), o.id.end());
    out.side.insert(out.side.end(), o.side.begin(), o.side.end());
    for (auto [dst, src] : {std::pair{&out.x, &o.x}, {&out.y, &o.y}, {&out.z, &o.z}, {&out.ux, &o.ux},
                            {&out.uy, &o.uy}, {&out.uz, &o.uz}, {&out.w, &o.w}, {&out.path, &o.path}})
      dst->insert(dst->end(), src->begin(), src->end());
    out.scatterings.insert(out.scatterings.end(), o.scatterings.begin(), o.scatterings.end());
  }
  return out;
}

}  // namespace scattuccino
