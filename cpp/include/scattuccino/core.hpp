// C++ core of Scattuccino (stage S5): photon transport in voxel geometry, mirroring
// src/scattuccino/voxels.py with the same counter-based random numbers (src/scattuccino/random.py)
// and the same order of draws per packet, so that CPU (and later GPU) results agree with the
// Python reference packet by packet up to the rounding of the elementary functions.
#pragma once

#include <cstdint>
#include <vector>

namespace scattuccino {

// Counter-based uniform numbers in (0, 1): u = mix(key(seed, id) ^ (k·G + 3)), 53 bits.
struct CounterRandom {
  std::uint64_t key;
  std::uint64_t counter = 0;
  CounterRandom(std::uint64_t seed, std::uint64_t id);
  double uniform();
};

enum class PhaseKind : int { HenyeyGreenstein = 0, ModifiedHG = 1, Gegenbauer = 2, Tabulated = 3 };

struct Phase {
  PhaseKind kind = PhaseKind::HenyeyGreenstein;
  double p1 = 0.0;  // g (HG), g_hg (MHG), g_gk (GK)
  double p2 = 0.0;  // beta (MHG), alpha (GK)
  std::vector<double> mu, p, cdf;  // tabulated (normalised as in Python)
  double sample(double xi) const;
};

struct Material {
  double mua = 0.0, mus = 0.0, n = 1.0;
  Phase phase;
};

struct VoxelInput {
  int nx = 1, ny = 1, nz = 1;
  double dx = 1.0, dy = 1.0, dz = 1.0;
  std::vector<std::int32_t> labels;  // (nx, ny, nz), C order
  std::vector<Material> materials;
  double n_above = 1.0, n_below = 1.0;
  bool periodic = true;
  double cos_incident = 1.0;
  std::uint64_t seed = 1;
  std::int64_t packets = 0;
  double w_threshold = 1e-4, roulette_m = 10.0;
  bool record = true;  // protocol of the leaving packets
  bool score_voxels = true;  // absorbed power per voxel (the GPU skips the atomics if false)
  // 0 double (reference); 1 single, position in float (A); 2 single, position in double (B)
  int precision = 0;
  int threads = 1;
};

struct VoxelOutput {
  double r_specular = 0.0;
  double sum_r = 0.0, sum_r2 = 0.0, sum_t = 0.0, sum_t2 = 0.0, sum_side = 0.0;  // exit weights
  std::int64_t n_exit = 0, n_roulette = 0;  // packets that left / were killed by the roulette
  std::vector<double> absorbed;         // per material, summed weight (not divided by packets)
  std::vector<double> absorbed_voxels;  // per voxel
  // protocol of leaving packets, in packet order
  std::vector<std::int64_t> id;
  std::vector<int> side;
  std::vector<double> x, y, ux, uy, uz, w, zmax;
  std::vector<double> path;  // (exits, materials)
  std::vector<std::int64_t> scatterings;
};

// Derived launch quantities (entry voxel, refracted direction, specular reflection).
struct Setup {
  double lo[2] = {0.0, 0.0};
  int i0 = 0, j0 = 0;
  double r_specular = 0.0, sin_t0 = 0.0, cos_t0 = 1.0;
};
Setup prepare(const VoxelInput& in);

VoxelOutput simulate_voxels(const VoxelInput& in);

// CUDA backend (stage S5b), available if built with SCATTUCCINO_ENABLE_CUDA.
int cuda_device_count();
VoxelOutput simulate_voxels_cuda(const VoxelInput& in, int device);

}  // namespace scattuccino
