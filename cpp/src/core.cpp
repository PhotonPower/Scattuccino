// Voxel transport, see cpp/include/scattuccino/core.hpp. Line-by-line mirror of voxels.py.
#include "scattuccino/core.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace scattuccino {

namespace {
constexpr std::uint64_t M1 = 0xBF58476D1CE4E5B9ULL, M2 = 0x94D049BB133111EBULL,
                        G = 0x9E3779B97F4A7C15ULL;
constexpr double PI = 3.141592653589793;

std::uint64_t mix(std::uint64_t x) {
  x = (x ^ (x >> 30)) * M1;
  x = (x ^ (x >> 27)) * M2;
  return x ^ (x >> 31);
}

// unpolarised Fresnel reflectance and cos θ_t (TIR: R = 1), as layers.fresnel_unpolarised
void fresnel(double n1, double n2, double cos_i, double& r, double& cos_t) {
  cos_i = std::abs(cos_i);
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
  cos_t = std::sqrt(std::clamp(1 - sin_t2, 0.0, 1.0));
  const double rs = (n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t);
  const double rp = (n2 * cos_i - n1 * cos_t) / (n2 * cos_i + n1 * cos_t);
  r = 0.5 * (rs * rs + rp * rp);
}

double hg(double g, double xi) {
  if (std::abs(g) < 1e-12) return 2 * xi - 1;
  const double t = (1 - g * g) / (1 - g + 2 * g * xi);
  return std::clamp((1 + g * g - t * t) / (2 * g), -1.0, 1.0);
}

// new direction after deflection, as layers._scatter
void scatter(double& ux, double& uy, double& uz, double cos_t, double phi) {
  const double sin_t = std::sqrt(std::clamp(1 - cos_t * cos_t, 0.0, 1.0));
  const double cp = std::cos(phi), sp = std::sin(phi);
  if (std::abs(uz) > 1 - 1e-12) {
    const double s = uz > 0 ? 1.0 : (uz < 0 ? -1.0 : 0.0);
    ux = sin_t * cp;
    uy = sin_t * sp;
    uz = s * cos_t;
    return;
  }
  const double den = std::sqrt(std::max(1 - uz * uz, 1e-300));
  const double nx = sin_t * (ux * uz * cp - uy * sp) / den + ux * cos_t;
  const double ny = sin_t * (uy * uz * cp + ux * sp) / den + uy * cos_t;
  const double nz = -sin_t * cp * den + uz * cos_t;
  ux = nx;
  uy = ny;
  uz = nz;
}
}  // namespace

CounterRandom::CounterRandom(std::uint64_t seed, std::uint64_t id)
    : key(mix(seed * G + 1) ^ mix(id * G + 2)) {}

double CounterRandom::uniform() {
  const std::uint64_t x = mix(key ^ (counter * G + 3));
  ++counter;
  return (static_cast<double>(x >> 11) + 0.5) * 0x1.0p-53;
}

double Phase::sample(double xi) const {
  switch (kind) {
    case PhaseKind::HenyeyGreenstein:
      return hg(p1, xi);
    case PhaseKind::ModifiedHG: {
      if (xi < p2) return hg(p1, xi / std::max(p2, 1e-300));
      return std::cbrt(2 * ((xi - p2) / std::max(1 - p2, 1e-300)) - 1);
    }
    case PhaseKind::Gegenbauer: {
      const double g = p1, a = p2;
      const double norm = std::pow(1 - g, -2 * a) - std::pow(1 + g, -2 * a);
      const double u = std::pow(std::pow(1 + g, -2 * a) + xi * norm, -1 / a);
      return std::clamp((1 + g * g - u) / (2 * g), -1.0, 1.0);
    }
    case PhaseKind::Tabulated: {
      const auto it = std::upper_bound(cdf.begin(), cdf.end(), xi);
      std::size_t k = static_cast<std::size_t>(std::max<std::ptrdiff_t>(it - cdf.begin() - 1, 0));
      k = std::min(k, mu.size() - 2);
      const double m0 = mu[k], m1 = mu[k + 1], q0 = p[k], q1 = p[k + 1];
      const double slope = (q1 - q0) / (m1 - m0);
      const double c = xi - cdf[k];
      double dq;
      if (std::abs(slope) < 1e-12 * std::max(q0, 1e-300))
        dq = c / (2 * PI * q0);
      else
        dq = (-q0 + std::sqrt(std::max(q0 * q0 + 2 * slope * c / (2 * PI), 0.0))) / slope;
      return std::clamp(m0 + dq, -1.0, 1.0);
    }
  }
  return 0.0;
}

VoxelOutput simulate_voxels(const VoxelInput& in) {
  const int nx = in.nx, ny = in.ny, nz = in.nz;
  const std::size_t nvox = static_cast<std::size_t>(nx) * ny * nz;
  if (in.labels.size() != nvox) throw std::invalid_argument("labels must have nx·ny·nz entries");
  const std::size_t nm = in.materials.size();
  for (auto l : in.labels)
    if (l < 0 || static_cast<std::size_t>(l) >= nm) throw std::invalid_argument("label out of range");
  const std::array<double, 3> size{in.dx, in.dy, in.dz};
  const std::array<double, 3> lo{-nx * in.dx / 2, -ny * in.dy / 2, 0.0};
  const std::array<int, 3> dims{nx, ny, nz};
  auto label = [&](int i, int j, int k) {
    return in.labels[(static_cast<std::size_t>(i) * ny + j) * nz + k];
  };
  const int i0 = std::min(static_cast<int>(-lo[0] / size[0]), nx - 1);
  const int j0 = std::min(static_cast<int>(-lo[1] / size[1]), ny - 1);
  const Material& entry = in.materials[static_cast<std::size_t>(label(i0, j0, 0))];
  double r_sp, cos_t0;
  fresnel(in.n_above, entry.n, in.cos_incident, r_sp, cos_t0);
  const double sin_t0 = in.n_above / entry.n * std::sqrt(1 - in.cos_incident * in.cos_incident);

  VoxelOutput out;
  out.r_specular = r_sp;
  out.absorbed.assign(nm, 0.0);
  out.absorbed_voxels.assign(nvox, 0.0);
  const std::int64_t n_packets = in.packets;
  const int threads = std::max(in.threads, 1);
  std::vector<VoxelOutput> parts(static_cast<std::size_t>(threads));
#pragma omp parallel for num_threads(threads) schedule(static)
  for (int t = 0; t < threads; ++t) {
    VoxelOutput& o = parts[static_cast<std::size_t>(t)];
    o.absorbed.assign(nm, 0.0);
    o.absorbed_voxels.assign(nvox, 0.0);
    const std::int64_t begin = n_packets * t / threads, end = n_packets * (t + 1) / threads;
    std::vector<double> path(nm);
    for (std::int64_t id = begin; id < end; ++id) {
      CounterRandom rng(in.seed, static_cast<std::uint64_t>(id));
      std::array<double, 3> pos{0.0, 0.0, 0.0}, u{sin_t0, 0.0, cos_t0};
      std::array<int, 3> idx{i0, j0, 0};
      double w = 1.0 - r_sp, s_left = 0.0, zmax = 0.0;
      std::int64_t nsc = 0;
      std::fill(path.begin(), path.end(), 0.0);
      while (true) {
        if (s_left <= 0) s_left = -std::log(rng.uniform());
        const int mat = label(idx[0], idx[1], idx[2]);
        const Material& M = in.materials[static_cast<std::size_t>(mat)];
        const double mt = M.mua + M.mus;
        const double step = mt > 0 ? s_left / mt : HUGE_VAL;
        std::array<double, 3> d_ax;
        for (int a = 0; a < 3; ++a) {
          const double f = lo[static_cast<std::size_t>(a)] + idx[static_cast<std::size_t>(a)] * size[static_cast<std::size_t>(a)];
          const double ua = u[static_cast<std::size_t>(a)];
          double d = ua > 0 ? (f + size[static_cast<std::size_t>(a)] - pos[static_cast<std::size_t>(a)]) / ua
                    : (ua < 0 ? (f - pos[static_cast<std::size_t>(a)]) / ua : HUGE_VAL);
          d_ax[static_cast<std::size_t>(a)] = std::max(d, 0.0);
        }
        int axis = 0;
        for (int a = 1; a < 3; ++a)
          if (d_ax[static_cast<std::size_t>(a)] < d_ax[static_cast<std::size_t>(axis)]) axis = a;
        const double d_face = d_ax[static_cast<std::size_t>(axis)];
        const bool hit = d_face <= step;
        const double dist = hit ? d_face : step;
        for (int a = 0; a < 3; ++a) pos[static_cast<std::size_t>(a)] += u[static_cast<std::size_t>(a)] * dist;
        path[static_cast<std::size_t>(mat)] += dist;
        zmax = std::max(zmax, pos[2]);
        s_left = hit ? s_left - d_face * mt : 0.0;
        if (hit) {
          const auto ax = static_cast<std::size_t>(axis);
          const int direction = u[ax] > 0 ? 1 : (u[ax] < 0 ? -1 : 0);
          std::array<int, 3> nb = idx;
          nb[ax] += direction;
          const double n1 = M.n;
          bool out_top = axis == 2 && nb[2] < 0, out_bot = axis == 2 && nb[2] >= nz;
          bool out_lat = axis < 2 && (nb[ax] < 0 || nb[ax] >= dims[ax]);
          if (in.periodic && out_lat) {
            nb[ax] = ((nb[ax] % dims[ax]) + dims[ax]) % dims[ax];
            pos[ax] -= direction * dims[ax] * size[ax];
            out_lat = false;
          }
          const bool inside = !(out_top || out_bot || out_lat);
          const double n2 = out_top ? in.n_above
                          : out_bot ? in.n_below
                          : out_lat ? n1
                                    : in.materials[static_cast<std::size_t>(label(nb[0], nb[1], nb[2]))].n;
          const double cos_i = std::abs(u[ax]);
          double r = 0.0, cos_t = cos_i;
          bool refl = false;
          if (n1 != n2) {
            fresnel(n1, n2, cos_i, r, cos_t);
            refl = rng.uniform() < r;
          }
          if (refl) {
            u[ax] = -u[ax];
          } else {
            const double ratio = n1 / n2;
            for (int a = 0; a < 3; ++a) u[static_cast<std::size_t>(a)] *= ratio;
            u[ax] = direction * cos_t;
            if (!inside) {
              if (in.record) {
                o.id.push_back(id);
                o.side.push_back(out_top ? -1 : (out_bot ? 1 : 0));
                o.x.push_back(pos[0]);
                o.y.push_back(pos[1]);
                o.ux.push_back(u[0]);
                o.uy.push_back(u[1]);
                o.uz.push_back(u[2]);
                o.w.push_back(w);
                o.zmax.push_back(zmax);
                o.path.insert(o.path.end(), path.begin(), path.end());
                o.scatterings.push_back(nsc);
              }
              break;
            }
            idx = nb;
          }
        } else {
          const double dw = mt > 0 ? w * (M.mua / mt) : 0.0;
          o.absorbed[static_cast<std::size_t>(mat)] += dw;
          o.absorbed_voxels[(static_cast<std::size_t>(idx[0]) * ny + idx[1]) * nz + idx[2]] += dw;
          w -= dw;
          const double ct = M.phase.sample(rng.uniform());
          const double phi = 2 * PI * rng.uniform();
          scatter(u[0], u[1], u[2], ct, phi);
          ++nsc;
          if (w < in.w_threshold) {
            if (rng.uniform() < 1.0 / in.roulette_m) {
              w *= in.roulette_m;
            } else {
              break;
            }
          }
        }
      }
    }
  }
  for (const auto& o : parts) {
    for (std::size_t k = 0; k < nm; ++k) out.absorbed[k] += o.absorbed[k];
    for (std::size_t k = 0; k < nvox; ++k) out.absorbed_voxels[k] += o.absorbed_voxels[k];
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
