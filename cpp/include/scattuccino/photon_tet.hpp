// The life of one photon packet in a tetrahedral mesh (stage S4b), CPU and (later) GPU.
//
// Same physics and order of random draws as the voxel walk (photon.hpp): the packet moves to the
// next interaction or the next tetrahedron face; at faces between different refractive indices
// (or at the mesh boundary, with a per-face exterior index) it is reflected or refracted by
// Fresnel about the face normal; faces between equal indices are crossed without drawing a random
// number, so refining a mesh of one material does not change the random stream.
#pragma once

#include "scattuccino/photon.hpp"

namespace scattuccino::photon {

struct TetMesh {
  const double* normal;     // (tets, 4, 3) outward unit normals of the faces
  const double* offset;     // (tets, 4) plane offsets: n·x = offset on the face
  const int* neighbour;     // (tets, 4) tetrahedron across the face, −1 at the boundary
  const double* exterior_n; // (tets, 4) index outside boundary faces
  const int* label;         // (tets) material
  const MaterialView<double>* materials;
  double w_threshold, roulette_m;
};

struct TetExit {
  int side;  // −1 leaving with u_z < 0, +1 with u_z ≥ 0, 2 roulette
  double x, y, z, ux, uy, uz, w;
  long long scatterings;
};

// Refraction of u through a face with outward normal n (cos_i = n·u > 0), indices n1 → n2,
// cos_t from Fresnel: u_t = (n1/n2) u + (cos_t − (n1/n2) cos_i) n.
SCAT_HD inline void refract(double* u, const double* n, double cos_i, double cos_t, double ratio) {
  for (int a = 0; a < 3; ++a) u[a] = ratio * u[a] + (cos_t - ratio * cos_i) * n[a];
}

template <class Deposit, class Path>
SCAT_HD TetExit trace_tet(const TetMesh& m, std::uint64_t seed, std::uint64_t id, const double* p0,
                          const double* u0, double w0, int tet0, Deposit deposit, Path add_path) {
  Rng rng(seed, id);
  double pos[3] = {p0[0], p0[1], p0[2]}, u[3] = {u0[0], u0[1], u0[2]};
  int tet = tet0;
  double w = w0, s_left = 0.0;
  long long nsc = 0;
  while (true) {
    if (s_left <= 0) s_left = -log(rng.uniform());
    const int mat = m.label[tet];
    const MaterialView<double>& M = m.materials[mat];
    const double mt = M.mua + M.mus;
    const double step = mt > 0 ? s_left / mt : HUGE_VAL;
    // exit face: smallest non-negative distance among faces the direction points out of
    int face = -1;
    double d_face = HUGE_VAL;
    for (int f = 0; f < 4; ++f) {
      const double* n = m.normal + (static_cast<long long>(tet) * 4 + f) * 3;
      const double un = n[0] * u[0] + n[1] * u[1] + n[2] * u[2];
      if (un <= 0) continue;
      const double d = dmax((m.offset[tet * 4 + f] - (n[0] * pos[0] + n[1] * pos[1] + n[2] * pos[2])) / un, 0.0);
      if (d < d_face) {
        d_face = d;
        face = f;
      }
    }
    const bool hit = face >= 0 && d_face <= step;
    const double dist = hit ? d_face : step;
    for (int a = 0; a < 3; ++a) pos[a] += u[a] * dist;
    add_path(mat, dist);
    s_left = hit ? s_left - d_face * mt : 0.0;
    if (hit) {
      const double* n = m.normal + (static_cast<long long>(tet) * 4 + face) * 3;
      const int nb = m.neighbour[tet * 4 + face];
      const double n1 = M.n;
      const double n2 = nb < 0 ? m.exterior_n[tet * 4 + face] : m.materials[m.label[nb]].n;
      const double cos_i = n[0] * u[0] + n[1] * u[1] + n[2] * u[2];
      double r = 0.0, cos_t = cos_i;
      bool refl = false;
      if (n1 != n2) {
        fresnel<double>(n1, n2, cos_i, r, cos_t);
        refl = rng.uniform() < r;
      }
      if (refl) {
        for (int a = 0; a < 3; ++a) u[a] -= 2 * cos_i * n[a];
      } else {
        if (n1 != n2) refract(u, n, cos_i, cos_t, n1 / n2);
        if (nb < 0) return TetExit{u[2] < 0 ? -1 : 1, pos[0], pos[1], pos[2], u[0], u[1], u[2], w, nsc};
        tet = nb;
      }
    } else {
      const double dw = mt > 0 ? w * (M.mua / mt) : 0.0;
      deposit(mat, tet, dw);
      w -= dw;
      const double ct = sample<double>(M.phase, rng.uniform());
      const double phi = 2 * PI * rng.uniform();
      scatter<double>(u[0], u[1], u[2], ct, phi);
      ++nsc;
      if (w < m.w_threshold) {
        if (rng.uniform() < 1.0 / m.roulette_m) {
          w *= m.roulette_m;
        } else {
          return TetExit{2, pos[0], pos[1], pos[2], u[0], u[1], u[2], 0.0, nsc};
        }
      }
    }
  }
}

}  // namespace scattuccino::photon
