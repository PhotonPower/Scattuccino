// Transport in tetrahedral meshes (stage S4b): host interface.
#pragma once

#include <cstdint>
#include <vector>

#include "scattuccino/core.hpp"

namespace scattuccino {

struct TetInput {
  std::vector<double> nodes;     // (nodes, 3) in mm
  std::vector<std::int32_t> tets;   // (tets, 4) node indices
  std::vector<std::int32_t> labels; // (tets) material
  std::vector<Material> materials;
  // exterior index for each boundary face, looked up by the face's outward normal: faces with
  // n_z < −0.5 get n_top, n_z > 0.5 n_bottom, others n_side (simple and enough for boxes and
  // spheres); or an explicit value for all boundary faces if n_all > 0.
  double n_top = 1.0, n_bottom = 1.0, n_side = 1.0, n_all = 0.0;
  std::vector<double> entry_point{0.0, 0.0, 0.0};  // on the boundary
  std::vector<double> direction{0.0, 0.0, 1.0};    // incident direction (outside the mesh)
  std::uint64_t seed = 1;
  std::int64_t packets = 0;
  double w_threshold = 1e-4, roulette_m = 10.0;
  bool record = true;
  int threads = 1;
};

struct TetOutput {
  double r_specular = 0.0;
  int entry_tet = -1;
  std::vector<double> absorbed;      // per material
  std::vector<double> absorbed_tets; // per tetrahedron
  std::vector<std::int64_t> id;
  std::vector<int> side;
  std::vector<double> x, y, z, ux, uy, uz, w, path;
  std::vector<std::int64_t> scatterings;
};

TetOutput simulate_tets(const TetInput& in);

}  // namespace scattuccino
