// The life of one photon packet in voxel geometry, shared by the CPU core (core.cpp) and the
// CUDA kernel (core_cuda.cu) so that both run literally the same code (stage S5).
//
// Mirrors src/scattuccino/voxels.py: counter-based random numbers per packet, the same order of
// draws, implicit capture, Russian roulette, Fresnel at voxel faces between different indices.
#pragma once

#include <cmath>
#include <cstdint>

#if defined(__CUDACC__)
#define SCAT_HD __host__ __device__
#else
#define SCAT_HD
#endif

namespace scattuccino::photon {

constexpr std::uint64_t M1 = 0xBF58476D1CE4E5B9ULL, M2 = 0x94D049BB133111EBULL,
                        G = 0x9E3779B97F4A7C15ULL;
constexpr double PI = 3.141592653589793;

SCAT_HD inline double dmin(double a, double b) { return a < b ? a : b; }
SCAT_HD inline double dmax(double a, double b) { return a > b ? a : b; }
SCAT_HD inline double clampd(double v, double lo, double hi) { return v < lo ? lo : (v > hi ? hi : v); }

SCAT_HD inline std::uint64_t mix(std::uint64_t x) {
  x = (x ^ (x >> 30)) * M1;
  x = (x ^ (x >> 27)) * M2;
  return x ^ (x >> 31);
}

struct Rng {
  std::uint64_t key, counter;
  SCAT_HD Rng(std::uint64_t seed, std::uint64_t id) : key(mix(seed * G + 1) ^ mix(id * G + 2)), counter(0) {}
  SCAT_HD double uniform() {
    const std::uint64_t x = mix(key ^ (counter * G + 3));
    ++counter;
    return (static_cast<double>(x >> 11) + 0.5) * 0x1.0p-53;
  }
};

// Phase function as plain data (tabulated arrays as pointers: host or device memory).
struct PhaseView {
  int kind = 0;  // 0 HG, 1 MHG, 2 Gegenbauer, 3 tabulated
  double p1 = 0.0, p2 = 0.0;
  const double* mu = nullptr;
  const double* p = nullptr;
  const double* cdf = nullptr;
  int n = 0;
};

struct MaterialView {
  double mua, mus, n;
  PhaseView phase;
};

SCAT_HD inline double hg(double g, double xi) {
  if (fabs(g) < 1e-12) return 2 * xi - 1;
  const double t = (1 - g * g) / (1 - g + 2 * g * xi);
  return clampd((1 + g * g - t * t) / (2 * g), -1.0, 1.0);
}

SCAT_HD inline double sample(const PhaseView& ph, double xi) {
  switch (ph.kind) {
    case 0:
      return hg(ph.p1, xi);
    case 1:
      if (xi < ph.p2) return hg(ph.p1, xi / dmax(ph.p2, 1e-300));
      return cbrt(2 * ((xi - ph.p2) / dmax(1 - ph.p2, 1e-300)) - 1);
    case 2: {
      const double g = ph.p1, a = ph.p2;
      const double norm = pow(1 - g, -2 * a) - pow(1 + g, -2 * a);
      const double u = pow(pow(1 + g, -2 * a) + xi * norm, -1 / a);
      return clampd((1 + g * g - u) / (2 * g), -1.0, 1.0);
    }
    default: {
      // first index with cdf > xi, minus one (as searchsorted(side="right") − 1)
      int lo = 0, hi = ph.n;
      while (lo < hi) {
        const int mid = (lo + hi) / 2;
        if (ph.cdf[mid] > xi) hi = mid; else lo = mid + 1;
      }
      int k = lo - 1;
      if (k < 0) k = 0;
      if (k > ph.n - 2) k = ph.n - 2;
      const double m0 = ph.mu[k], m1 = ph.mu[k + 1], q0 = ph.p[k], q1 = ph.p[k + 1];
      const double slope = (q1 - q0) / (m1 - m0);
      const double c = xi - ph.cdf[k];
      double dq;
      if (fabs(slope) < 1e-12 * dmax(q0, 1e-300))
        dq = c / (2 * PI * q0);
      else
        dq = (-q0 + sqrt(dmax(q0 * q0 + 2 * slope * c / (2 * PI), 0.0))) / slope;
      return clampd(m0 + dq, -1.0, 1.0);
    }
  }
}

SCAT_HD inline void fresnel(double n1, double n2, double cos_i, double& r, double& cos_t) {
  cos_i = fabs(cos_i);
  if (n1 == n2) {
    r = 0.0;
    cos_t = cos_i;
    return;
  }
  const double sin_t2 = (n1 / n2) * (n1 / n2) * (1 - cos_i * cos_i);
  if (sin_t2 >= 1.0) {
    r = 1.0;
    cos_t = 0.0;
    return;
  }
  cos_t = sqrt(clampd(1 - sin_t2, 0.0, 1.0));
  const double rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t);
  const double rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t);
  r = 0.5 * (rs * rs + rp * rp);
}

SCAT_HD inline void scatter(double& ux, double& uy, double& uz, double cos_t, double phi) {
  const double sin_t = sqrt(clampd(1 - cos_t * cos_t, 0.0, 1.0));
  const double cp = cos(phi), sp = sin(phi);
  if (fabs(uz) > 1 - 1e-12) {
    const double s = uz > 0 ? 1.0 : (uz < 0 ? -1.0 : 0.0);
    ux = sin_t * cp;
    uy = sin_t * sp;
    uz = s * cos_t;
    return;
  }
  const double den = sqrt(dmax(1 - uz * uz, 1e-300));
  const double nx = sin_t * (ux * uz * cp - uy * sp) / den + ux * cos_t;
  const double ny = sin_t * (uy * uz * cp + ux * sp) / den + uy * cos_t;
  const double nz = -sin_t * cp * den + uz * cos_t;
  ux = nx;
  uy = ny;
  uz = nz;
}

struct Geometry {
  int nx, ny, nz;
  double size[3], lo[3];
  const std::int32_t* labels;
  const MaterialView* materials;
  double n_above, n_below;
  bool periodic;
  double w_threshold, roulette_m;
  SCAT_HD int label(int i, int j, int k) const {
    return labels[(static_cast<long long>(i) * ny + j) * nz + k];
  }
};

// Exit of a packet.
struct Exit {
  int side;  // −1 top, +1 bottom, 0 side, 2 killed by roulette (not recorded as exit)
  double x, y, ux, uy, uz, w, zmax;
  long long scatterings;
};

// Runs one packet; Deposit(material, i, j, k, dw) and Path(material, dist) are callables.
template <class Deposit, class Path>
SCAT_HD Exit trace(const Geometry& g, std::uint64_t seed, std::uint64_t id, double sin_t0,
                   double cos_t0, double w0, int i0, int j0, Deposit deposit, Path add_path) {
  Rng rng(seed, id);
  double pos[3] = {0.0, 0.0, 0.0}, u[3] = {sin_t0, 0.0, cos_t0};
  int idx[3] = {i0, j0, 0};
  const int dims[3] = {g.nx, g.ny, g.nz};
  double w = w0, s_left = 0.0, zmax = 0.0;
  long long nsc = 0;
  while (true) {
    if (s_left <= 0) s_left = -log(rng.uniform());
    const int mat = g.label(idx[0], idx[1], idx[2]);
    const MaterialView& M = g.materials[mat];
    const double mt = M.mua + M.mus;
    const double step = mt > 0 ? s_left / mt : HUGE_VAL;
    double d_ax[3];
    for (int a = 0; a < 3; ++a) {
      const double f = g.lo[a] + idx[a] * g.size[a];
      const double ua = u[a];
      const double d = ua > 0 ? (f + g.size[a] - pos[a]) / ua : (ua < 0 ? (f - pos[a]) / ua : HUGE_VAL);
      d_ax[a] = dmax(d, 0.0);
    }
    int axis = 0;
    if (d_ax[1] < d_ax[axis]) axis = 1;
    if (d_ax[2] < d_ax[axis]) axis = 2;
    const double d_face = d_ax[axis];
    const bool hit = d_face <= step;
    const double dist = hit ? d_face : step;
    for (int a = 0; a < 3; ++a) pos[a] += u[a] * dist;
    add_path(mat, dist);
    zmax = dmax(zmax, pos[2]);
    s_left = hit ? s_left - d_face * mt : 0.0;
    if (hit) {
      const int direction = u[axis] > 0 ? 1 : (u[axis] < 0 ? -1 : 0);
      int nb[3] = {idx[0], idx[1], idx[2]};
      nb[axis] += direction;
      const double n1 = M.n;
      const bool out_top = axis == 2 && nb[2] < 0, out_bot = axis == 2 && nb[2] >= g.nz;
      bool out_lat = axis < 2 && (nb[axis] < 0 || nb[axis] >= dims[axis]);
      if (g.periodic && out_lat) {
        nb[axis] = ((nb[axis] % dims[axis]) + dims[axis]) % dims[axis];
        pos[axis] -= direction * dims[axis] * g.size[axis];
        out_lat = false;
      }
      const bool inside = !(out_top || out_bot || out_lat);
      const double n2 = out_top ? g.n_above
                      : out_bot ? g.n_below
                      : out_lat ? n1
                                : g.materials[g.label(nb[0], nb[1], nb[2])].n;
      const double cos_i = fabs(u[axis]);
      double r = 0.0, cos_t = cos_i;
      bool refl = false;
      if (n1 != n2) {
        fresnel(n1, n2, cos_i, r, cos_t);
        refl = rng.uniform() < r;
      }
      if (refl) {
        u[axis] = -u[axis];
      } else {
        const double ratio = n1 / n2;
        for (int a = 0; a < 3; ++a) u[a] *= ratio;
        u[axis] = direction * cos_t;
        if (!inside) {
          return Exit{out_top ? -1 : (out_bot ? 1 : 0), pos[0], pos[1], u[0], u[1], u[2], w, zmax, nsc};
        }
        idx[0] = nb[0];
        idx[1] = nb[1];
        idx[2] = nb[2];
      }
    } else {
      const double dw = mt > 0 ? w * (M.mua / mt) : 0.0;
      deposit(mat, idx[0], idx[1], idx[2], dw);
      w -= dw;
      const double ct = sample(M.phase, rng.uniform());
      const double phi = 2 * PI * rng.uniform();
      scatter(u[0], u[1], u[2], ct, phi);
      ++nsc;
      if (w < g.w_threshold) {
        if (rng.uniform() < 1.0 / g.roulette_m) {
          w *= g.roulette_m;
        } else {
          return Exit{2, pos[0], pos[1], u[0], u[1], u[2], 0.0, zmax, nsc};
        }
      }
    }
  }
}

}  // namespace scattuccino::photon
