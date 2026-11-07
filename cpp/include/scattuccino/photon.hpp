// The life of one photon packet in voxel geometry, shared by the CPU core (core.cpp) and the
// CUDA kernel (core_cuda.cu) so that both run literally the same code (stage S5).
//
// Mirrors src/scattuccino/voxels.py: counter-based random numbers per packet, the same order of
// draws, implicit capture, Russian roulette, Fresnel at voxel faces between different indices.
#pragma once

#include <cmath>
#include <cstdint>
#include <type_traits>

#if defined(__CUDACC__)
#define SCAT_HD __host__ __device__
#else
#define SCAT_HD
#endif

namespace scattuccino::photon {

constexpr std::uint64_t M1 = 0xBF58476D1CE4E5B9ULL, M2 = 0x94D049BB133111EBULL,
                        G = 0x9E3779B97F4A7C15ULL;
constexpr double PI = 3.141592653589793;

// Precision of the packet code. Real (double or float): direction, step length, distances to the
// voxel faces, phase function, Fresnel, scattering. Always double: the packet weight, deposited
// energy, all sums, the protocol values and the random number itself (53 bit; converted to Real
// only when used, so the random stream of a packet is the same in both precisions). The position
// is Real too (measured against a double position with float direction: same statistics, 2.3 to
// 2.6 times slower on the GPU). With Real = double all functions are the former double code, bit
// for bit.
template <class R> constexpr bool is_single = std::is_same_v<R, float>;

SCAT_HD inline double dmin(double a, double b) { return a < b ? a : b; }
SCAT_HD inline double dmax(double a, double b) { return a > b ? a : b; }
template <class R> SCAT_HD inline R rmax(R a, R b) { return a > b ? a : b; }
template <class R> SCAT_HD inline R clampv(R v, R lo, R hi) { return v < lo ? lo : (v > hi ? hi : v); }

template <class R> SCAT_HD inline R rabs(R x) {
  if constexpr (is_single<R>) return fabsf(x); else return fabs(x);
}
template <class R> SCAT_HD inline R rsqrt_(R x) {
  if constexpr (is_single<R>) return sqrtf(x); else return sqrt(x);
}
template <class R> SCAT_HD inline R rlog(R x) {
  if constexpr (is_single<R>) return logf(x); else return log(x);
}
template <class R> SCAT_HD inline R rsin(R x) {
  if constexpr (is_single<R>) return sinf(x); else return sin(x);
}
template <class R> SCAT_HD inline R rcos(R x) {
  if constexpr (is_single<R>) return cosf(x); else return cos(x);
}
template <class R> SCAT_HD inline R rcbrt(R x) {
  if constexpr (is_single<R>) return cbrtf(x); else return cbrt(x);
}
template <class R> SCAT_HD inline R rpow(R x, R y) {
  if constexpr (is_single<R>) return powf(x, y); else return pow(x, y);
}
template <class R> SCAT_HD inline R infinity() {
  if constexpr (is_single<R>) return HUGE_VALF; else return HUGE_VAL;
}
// smallest guard values: below these g (HG) / slopes / polar distances count as zero
template <class R> SCAT_HD constexpr R tiny() { return is_single<R> ? R(1e-30) : R(1e-300); }
template <class R> SCAT_HD constexpr R small_g() { return is_single<R> ? R(1e-4) : R(1e-12); }
template <class R> SCAT_HD constexpr R pole() { return is_single<R> ? R(1e-7) : R(1e-12); }

// uniform number (double, 53 bit) as Real; in single precision never exactly 1 (log(1) = 0)
template <class R> SCAT_HD inline R as_real(double u) {
  if constexpr (is_single<R>) {
    const float f = static_cast<float>(u);
    return f < 1.0f ? f : 0x1.fffffep-1f;
  } else {
    return u;
  }
}

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
template <class R> struct PhaseView {
  int kind = 0;  // 0 HG, 1 MHG, 2 Gegenbauer, 3 tabulated
  R p1 = 0, p2 = 0;
  const R* mu = nullptr;
  const R* p = nullptr;
  const R* cdf = nullptr;
  int n = 0;
};

template <class R> struct MaterialView {
  R mua, mus, n;
  PhaseView<R> phase;
  double absfrac;  // mua / (mua + mus) in double: the weight deposited per interaction
};

template <class R> SCAT_HD inline R hg(R g, R xi) {
  if (rabs(g) < small_g<R>()) return 2 * xi - 1;
  const R t = (1 - g * g) / (1 - g + 2 * g * xi);
  return clampv((1 + g * g - t * t) / (2 * g), R(-1), R(1));
}

template <class R> SCAT_HD inline R sample(const PhaseView<R>& ph, R xi) {
  switch (ph.kind) {
    case 0:
      return hg(ph.p1, xi);
    case 1:
      if (xi < ph.p2) return hg(ph.p1, xi / rmax(ph.p2, tiny<R>()));
      return rcbrt(2 * ((xi - ph.p2) / rmax(1 - ph.p2, tiny<R>())) - 1);
    case 2: {
      const R g = ph.p1, a = ph.p2;
      const R norm = rpow(1 - g, -2 * a) - rpow(1 + g, -2 * a);
      const R u = rpow(rpow(1 + g, -2 * a) + xi * norm, -1 / a);
      return clampv((1 + g * g - u) / (2 * g), R(-1), R(1));
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
      const R m0 = ph.mu[k], m1 = ph.mu[k + 1], q0 = ph.p[k], q1 = ph.p[k + 1];
      const R slope = (q1 - q0) / (m1 - m0);
      const R c = xi - ph.cdf[k];
      R dq;
      if constexpr (is_single<R>) {
        // (−q0 + sqrt(q0² + 2·slope·c/2π)) / slope without the cancellation for small slopes
        const R s = rsqrt_(rmax(q0 * q0 + 2 * slope * c / R(2 * PI), R(0)));
        const R den = q0 + s;
        dq = den > tiny<R>() ? c / (R(PI) * den) : R(0);
      } else {
        if (rabs(slope) < small_g<R>() * rmax(q0, tiny<R>()))
          dq = c / (2 * PI * q0);
        else
          dq = (-q0 + rsqrt_(rmax(q0 * q0 + 2 * slope * c / (2 * PI), R(0)))) / slope;
      }
      return clampv(m0 + dq, R(-1), R(1));
    }
  }
}

template <class R> SCAT_HD inline void fresnel(R n1, R n2, R cos_i, R& r, R& cos_t) {
  cos_i = rabs(cos_i);
  if (n1 == n2) {
    r = 0;
    cos_t = cos_i;
    return;
  }
  const R sin_t2 = (n1 / n2) * (n1 / n2) * (1 - cos_i * cos_i);
  if (sin_t2 >= 1) {
    r = 1;
    cos_t = 0;
    return;
  }
  cos_t = rsqrt_(clampv(1 - sin_t2, R(0), R(1)));
  const R rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t);
  const R rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t);
  r = R(0.5) * (rs * rs + rp * rp);
}

template <class R> SCAT_HD inline void scatter(R& ux, R& uy, R& uz, R cos_t, R phi) {
  const R sin_t = rsqrt_(clampv(1 - cos_t * cos_t, R(0), R(1)));
  const R cp = rcos(phi), sp = rsin(phi);
  if (rabs(uz) > 1 - pole<R>()) {
    const R s = uz > 0 ? R(1) : (uz < 0 ? R(-1) : R(0));
    ux = sin_t * cp;
    uy = sin_t * sp;
    uz = s * cos_t;
    return;
  }
  R den2;
  if constexpr (is_single<R>) den2 = (1 - uz) * (1 + uz);  // exact for uz near ±1
  else den2 = 1 - uz * uz;
  const R den = rsqrt_(rmax(den2, tiny<R>()));
  R nx = sin_t * (ux * uz * cp - uy * sp) / den + ux * cos_t;
  R ny = sin_t * (uy * uz * cp + ux * sp) / den + uy * cos_t;
  R nz = -sin_t * cp * den + uz * cos_t;
  if constexpr (is_single<R>) {  // keep the direction a unit vector
    const R inv = 1 / rsqrt_(nx * nx + ny * ny + nz * nz);
    nx *= inv;
    ny *= inv;
    nz *= inv;
  }
  ux = nx;
  uy = ny;
  uz = nz;
}

template <class R> struct Geometry {
  int nx, ny, nz;
  R size[3], lo[3];
  const std::int32_t* labels;
  const MaterialView<R>* materials;
  R n_above, n_below;
  bool periodic;
  double w_threshold, roulette_m;
  SCAT_HD int label(int i, int j, int k) const {
    return labels[(static_cast<long long>(i) * ny + j) * nz + k];
  }
};

// Exit of a packet (always double).
struct Exit {
  int side;  // −1 top, +1 bottom, 0 side, 2 killed by roulette (not recorded as exit)
  double x, y, ux, uy, uz, w, zmax;
  long long scatterings;
};

// Runs one packet; Deposit(material, i, j, k, dw) and Path(material, dist) are callables.
template <class Real, class Deposit, class Path>
SCAT_HD Exit trace(const Geometry<Real>& g, std::uint64_t seed, std::uint64_t id, double sin_t0,
                   double cos_t0, double w0, int i0, int j0, Deposit deposit, Path add_path) {
  Rng rng(seed, id);
  Real pos[3] = {0, 0, 0};
  Real u[3] = {static_cast<Real>(sin_t0), 0, static_cast<Real>(cos_t0)};
  int idx[3] = {i0, j0, 0};
  const int dims[3] = {g.nx, g.ny, g.nz};
  double w = w0, zmax = 0.0;
  Real s_left = 0;
  long long nsc = 0;
  while (true) {
    if (s_left <= 0) s_left = -rlog(as_real<Real>(rng.uniform()));
    const int mat = g.label(idx[0], idx[1], idx[2]);
    const MaterialView<Real>& M = g.materials[mat];
    const Real mt = M.mua + M.mus;
    const Real step = mt > 0 ? s_left / mt : infinity<Real>();
    Real d_ax[3];
    for (int a = 0; a < 3; ++a) {
      const Real f = static_cast<Real>(g.lo[a]) + idx[a] * static_cast<Real>(g.size[a]);
      const Real ua = u[a];
      Real d;
      if (ua > 0) d = static_cast<Real>((f + static_cast<Real>(g.size[a])) - pos[a]) / ua;
      else if (ua < 0) d = static_cast<Real>(f - pos[a]) / ua;
      else d = infinity<Real>();
      d_ax[a] = rmax(d, Real(0));
    }
    int axis = 0;
    if (d_ax[1] < d_ax[axis]) axis = 1;
    if (d_ax[2] < d_ax[axis]) axis = 2;
    const Real d_face = d_ax[axis];
    const bool hit = d_face <= step;
    const Real dist = hit ? d_face : step;
    for (int a = 0; a < 3; ++a) pos[a] += static_cast<Real>(u[a]) * static_cast<Real>(dist);
    add_path(mat, dist);
    zmax = dmax(zmax, static_cast<double>(pos[2]));
    s_left = hit ? s_left - d_face * mt : Real(0);
    if (hit) {
      const int direction = u[axis] > 0 ? 1 : (u[axis] < 0 ? -1 : 0);
      if constexpr (is_single<Real>) {
        // the packet is on the face it crosses: remove the rounding of the step (the voxel
        // index then always belongs to the position)
        pos[axis] = static_cast<Real>(g.lo[axis]) +
                    (idx[axis] + (direction > 0 ? 1 : 0)) * static_cast<Real>(g.size[axis]);
      }
      int nb[3] = {idx[0], idx[1], idx[2]};
      nb[axis] += direction;
      const Real n1 = M.n;
      const bool out_top = axis == 2 && nb[2] < 0, out_bot = axis == 2 && nb[2] >= g.nz;
      bool out_lat = axis < 2 && (nb[axis] < 0 || nb[axis] >= dims[axis]);
      if (g.periodic && out_lat) {
        nb[axis] = ((nb[axis] % dims[axis]) + dims[axis]) % dims[axis];
        pos[axis] -= direction * dims[axis] * static_cast<Real>(g.size[axis]);
        out_lat = false;
      }
      const bool inside = !(out_top || out_bot || out_lat);
      const Real n2 = out_top ? g.n_above
                    : out_bot ? g.n_below
                    : out_lat ? n1
                              : g.materials[g.label(nb[0], nb[1], nb[2])].n;
      const Real cos_i = rabs(u[axis]);
      Real r = 0, cos_t = cos_i;
      bool refl = false;
      if (n1 != n2) {
        fresnel(n1, n2, cos_i, r, cos_t);
        refl = rng.uniform() < static_cast<double>(r);
      }
      if (refl) {
        u[axis] = -u[axis];
      } else {
        const Real ratio = n1 / n2;
        for (int a = 0; a < 3; ++a) u[a] *= ratio;
        u[axis] = direction * cos_t;
        if (!inside) {
          return Exit{out_top ? -1 : (out_bot ? 1 : 0),
                      static_cast<double>(pos[0]), static_cast<double>(pos[1]),
                      static_cast<double>(u[0]), static_cast<double>(u[1]), static_cast<double>(u[2]),
                      w, zmax, nsc};
        }
        idx[0] = nb[0];
        idx[1] = nb[1];
        idx[2] = nb[2];
      }
    } else {
      const double dw = mt > 0 ? w * M.absfrac : 0.0;
      deposit(mat, idx[0], idx[1], idx[2], dw);
      w -= dw;
      const Real ct = sample(M.phase, as_real<Real>(rng.uniform()));
      const Real phi = static_cast<Real>(2 * PI * rng.uniform());
      scatter(u[0], u[1], u[2], ct, phi);
      ++nsc;
      if (w < g.w_threshold) {
        if (rng.uniform() < 1.0 / g.roulette_m) {
          w *= g.roulette_m;
        } else {
          return Exit{2, static_cast<double>(pos[0]), static_cast<double>(pos[1]),
                      static_cast<double>(u[0]), static_cast<double>(u[1]), static_cast<double>(u[2]),
                      0.0, zmax, nsc};
        }
      }
    }
  }
}

}  // namespace scattuccino::photon
