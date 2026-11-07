// Polarised packets in voxel geometry (stage S6c, CPU): line-by-line mirror of
// src/scattuccino/polarisation.py (scatter_step, apply_scattering, interface, simulate_polarised)
// with the same counter-based random numbers and order of draws per packet.
//
// Host only for now; the GPU version follows with the FP32 option (D7).
#pragma once

#include <cmath>
#include <complex>

#include "scattuccino/photon.hpp"

namespace scattuccino::photon {

struct MuellerView {
  const double* mu = nullptr;
  const double* m12 = nullptr;
  const double* m33 = nullptr;
  const double* m34 = nullptr;
  int n = 0;
  bool pure = false;
};

struct PolMaterialView {
  double mua, mus, n;
  PhaseView<double> phase;  // tabulated S11 for θ
  MuellerView mueller;
};

// numpy.interp
inline double interp(double x, const double* xp, const double* fp, int n) {
  if (x < xp[0]) return fp[0];
  if (x > xp[n - 1]) return fp[n - 1];
  int lo = 0, hi = n;  // largest j with xp[j] <= x
  while (hi - lo > 1) {
    const int mid = (lo + hi) / 2;
    if (xp[mid] <= x) lo = mid; else hi = mid;
  }
  const int j = lo;
  if (j == n - 1) return fp[j];
  if (xp[j] == x) return fp[j];
  const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
  return slope * (x - xp[j]) + fp[j];
}

inline void ratios(const MuellerView& m, double mu, double& m12, double& m33, double& m34) {
  m12 = interp(mu, m.mu, m.m12, m.n);
  m33 = interp(mu, m.mu, m.m33, m.n);
  m34 = interp(mu, m.mu, m.m34, m.n);
  if (m.pure) {
    const double scale = std::sqrt(dmax(1 - m12 * m12, 0.0) / dmax(m33 * m33 + m34 * m34, 1e-300));
    m33 *= scale;
    m34 *= scale;
  }
}

inline void cross(const double* a, const double* b, double* c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
inline double dot3(const double* a, const double* b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline double norm3(const double* a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }

// apply_scattering (polarisation.py)
inline void apply_scattering(const MuellerView& mv, double* u, double* v, double* S, double mu, double phi) {
  const double q = S[1], uu = S[2];
  double m12, m33, m34;
  ratios(mv, mu, m12, m33, m34);
  double w[3];
  cross(u, v, w);
  const double cp = std::cos(phi), sp = std::sin(phi);
  double v1[3];
  for (int a = 0; a < 3; ++a) v1[a] = cp * v[a] + sp * w[a];
  const double c2 = std::cos(2 * phi), s2 = std::sin(2 * phi);
  const double q1 = q * c2 + uu * s2, u1 = -q * s2 + uu * c2;
  const double st = std::sqrt(clampv<double>(1 - mu * mu, 0.0, 1.0)), ct = mu;
  double un[3], vn[3];
  for (int a = 0; a < 3; ++a) {
    un[a] = ct * u[a] + st * v1[a];
    vn[a] = -st * u[a] + ct * v1[a];
  }
  const double i2 = 1.0 + m12 * q1, q2 = m12 + q1;
  const double u2 = m33 * u1 + m34 * S[3], v2 = -m34 * u1 + m33 * S[3];
  S[0] = 1.0;
  S[1] = q2 / i2;
  S[2] = u2 / i2;
  S[3] = v2 / i2;
  const double nu = norm3(un);
  for (double& x : un) x /= nu;
  const double d = dot3(vn, un);
  for (int a = 0; a < 3; ++a) vn[a] -= d * un[a];
  const double nv = norm3(vn);
  for (int a = 0; a < 3; ++a) {
    u[a] = un[a];
    v[a] = vn[a] / nv;
  }
}

// scatter_step: θ from S11, φ by rejection (pairs of draws), then apply_scattering
inline void scatter_pol(const PolMaterialView& M, double* u, double* v, double* S, Rng& rng) {
  const double mu = sample<double>(M.phase, rng.uniform());
  double m12, m33, m34;
  ratios(M.mueller, mu, m12, m33, m34);
  const double q = S[1], uu = S[2];
  const double bound = 1.0 + std::abs(m12) * std::sqrt(q * q + uu * uu);
  double phi;
  while (true) {
    const double trial = 2 * PI * rng.uniform();
    const double f = 1.0 + m12 * (q * std::cos(2 * trial) + uu * std::sin(2 * trial));
    if (rng.uniform() * bound <= f) {
      phi = trial;
      break;
    }
  }
  apply_scattering(M.mueller, u, v, S, mu, phi);
}

using cplx = std::complex<double>;

inline void fresnel_amplitudes(double n1, double n2, double cos_i, cplx& rp, cplx& rs, cplx& tp, cplx& ts, cplx& ct) {
  cos_i = std::abs(cos_i);
  const double sin2_t = (n1 / n2) * (n1 / n2) * (1 - cos_i * cos_i);
  ct = sin2_t <= 1 ? cplx(std::sqrt(dmax(1 - sin2_t, 0.0)), 0.0) : cplx(0.0, std::sqrt(dmax(sin2_t - 1, 0.0)));
  rs = (n1 * cos_i - n2 * ct) / (n1 * cos_i + n2 * ct);
  rp = (n2 * cos_i - n1 * ct) / (n2 * cos_i + n1 * ct);
  ts = 2 * n1 * cos_i / (n1 * cos_i + n2 * ct);
  tp = 2 * n1 * cos_i / (n2 * cos_i + n1 * ct);
}

inline double abs2(const cplx& z) {
  const double a = std::abs(z);
  return a * a;
}

inline void apply_amplitudes(const cplx& ap, const cplx& as, const double* in, double* out) {
  const double pp = abs2(ap), ss = abs2(as);
  const cplx c = ap * std::conj(as);
  out[0] = 0.5 * (pp + ss) * in[0] + 0.5 * (pp - ss) * in[1];
  out[1] = 0.5 * (pp - ss) * in[0] + 0.5 * (pp + ss) * in[1];
  out[2] = c.real() * in[2] + c.imag() * in[3];
  out[3] = -c.imag() * in[2] + c.real() * in[3];
}

// frame_to_plane: returns p̂, ŝ and the rotated Stokes vector
inline void frame_to_plane(const double* u, const double* v, const double* S, const double* normal, double* p_hat,
                           double* s_hat, double* Sr) {
  double w[3];
  cross(u, v, w);
  cross(u, normal, s_hat);
  const double nrm = norm3(s_hat);
  if (nrm < 1e-12) {
    for (int a = 0; a < 3; ++a) s_hat[a] = w[a];
  } else {
    for (int a = 0; a < 3; ++a) s_hat[a] /= nrm;
  }
  cross(s_hat, u, p_hat);
  const double phi = std::atan2(dot3(w, p_hat), dot3(v, p_hat));
  const double c2 = std::cos(2 * phi), s2 = std::sin(2 * phi);
  Sr[0] = S[0];
  Sr[1] = S[1] * c2 + S[2] * s2;
  Sr[2] = -S[1] * s2 + S[2] * c2;
  Sr[3] = S[3];
}

// interface (polarisation.py): returns true if reflected
inline bool interface_pol(double* u, double* v, double* S, double n1, double n2, const double* normal, double xi) {
  if (n1 == n2) return false;
  double p_hat[3], s_hat[3], st[4];
  frame_to_plane(u, v, S, normal, p_hat, s_hat, st);
  const double cos_i = std::abs(dot3(u, normal));
  cplx rp, rs, tp, ts, ct;
  fresnel_amplitudes(n1, n2, cos_i, rp, rs, tp, ts, ct);
  const double big_r = 0.5 * (abs2(rp) + abs2(rs)) + 0.5 * (abs2(rp) - abs2(rs)) * st[1];
  const bool refl = xi < big_r;
  double un[3];
  if (refl) {
    const double d = dot3(u, normal);
    for (int a = 0; a < 3; ++a) un[a] = u[a] - 2 * d * normal[a];
  } else {
    const double ratio = n1 / n2;
    const double c = ct.real();
    for (int a = 0; a < 3; ++a) un[a] = ratio * u[a] + (c - ratio * cos_i) * normal[a];
  }
  double sn[4];
  if (refl) apply_amplitudes(rp, rs, st, sn); else apply_amplitudes(tp, ts, st, sn);
  const double i0 = sn[0];
  for (int k = 0; k < 4; ++k) S[k] = sn[k] / i0;
  const double nu = norm3(un);
  for (int a = 0; a < 3; ++a) u[a] = un[a] / nu;
  cross(s_hat, u, v);
  return refl;
}

struct PolExit {
  int side;
  double x, y, u[3], v[3], S[4], w;
  long long scatterings;
};

struct PolGeometry {
  int nx, ny, nz;
  double size[3], lo[3];
  const std::int32_t* labels;
  const PolMaterialView* materials;
  double n_above, n_below;
  bool periodic;
  double w_threshold, roulette_m;
  int label(int i, int j, int k) const { return labels[(static_cast<long long>(i) * ny + j) * nz + k]; }
};

// One polarised packet; entry state (u0, v0, S0, w0) from the polarised entry (see core.cpp).
template <class Deposit>
PolExit trace_pol(const PolGeometry& g, std::uint64_t seed, std::uint64_t id, const double* u0, const double* v0,
                  const double* S0, double w0, int i0, int j0, Deposit deposit) {
  Rng rng(seed, id);
  double pos[3] = {0.0, 0.0, 0.0}, u[3] = {u0[0], u0[1], u0[2]}, v[3] = {v0[0], v0[1], v0[2]};
  double S[4] = {S0[0], S0[1], S0[2], S0[3]};
  int idx[3] = {i0, j0, 0};
  const int dims[3] = {g.nx, g.ny, g.nz};
  double w = w0, s_left = 0.0;
  long long nsc = 0;
  while (true) {
    if (s_left <= 0) s_left = -std::log(rng.uniform());
    const int mat = g.label(idx[0], idx[1], idx[2]);
    const PolMaterialView& M = g.materials[mat];
    const double mt = M.mua + M.mus;
    const double step = M.mus > 0 ? s_left / mt : HUGE_VAL;
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
    if (M.mus == 0 && M.mua > 0) w *= std::exp(-M.mua * dist);  // as the reference without scattering
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
      double normal[3] = {0.0, 0.0, 0.0};
      normal[axis] = direction;
      const bool refl = n1 != n2 && interface_pol(u, v, S, n1, n2, normal, rng.uniform());
      if (!refl) {
        if (!inside) {
          PolExit e{out_top ? -1 : (out_bot ? 1 : 0), pos[0], pos[1], {u[0], u[1], u[2]}, {v[0], v[1], v[2]},
                    {S[0], S[1], S[2], S[3]}, w, nsc};
          return e;
        }
        idx[0] = nb[0];
        idx[1] = nb[1];
        idx[2] = nb[2];
      }
    } else if (M.mus > 0) {
      const double factor = M.mus / mt;
      deposit(mat, idx[0], idx[1], idx[2], w * (1 - factor));
      w *= factor;
      scatter_pol(M, u, v, S, rng);
      ++nsc;
      if (w < g.w_threshold) {
        if (rng.uniform() < 1.0 / g.roulette_m) {
          w *= g.roulette_m;
        } else {
          PolExit e{2, pos[0], pos[1], {u[0], u[1], u[2]}, {v[0], v[1], v[2]}, {S[0], S[1], S[2], S[3]}, 0.0, nsc};
          return e;
        }
      }
    }
  }
}

}  // namespace scattuccino::photon
