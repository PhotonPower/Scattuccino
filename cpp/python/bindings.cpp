// Python module scattuccino._core (pybind11).
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scattuccino/core.hpp"

namespace py = pybind11;
using namespace scattuccino;

template <class T>
py::array_t<T> arr(const std::vector<T>& v) {
  py::array_t<T> a(static_cast<py::ssize_t>(v.size()));
  std::copy(v.begin(), v.end(), a.mutable_data());
  return a;
}

PYBIND11_MODULE(_core, m) {
  m.doc() = "C++ core of Scattuccino (S5): voxel transport, same random numbers as the Python reference.";
  py::enum_<PhaseKind>(m, "PhaseKind")
      .value("HenyeyGreenstein", PhaseKind::HenyeyGreenstein)
      .value("ModifiedHG", PhaseKind::ModifiedHG)
      .value("Gegenbauer", PhaseKind::Gegenbauer)
      .value("Tabulated", PhaseKind::Tabulated);
  py::class_<Phase>(m, "Phase")
      .def(py::init<>())
      .def_readwrite("kind", &Phase::kind)
      .def_readwrite("p1", &Phase::p1)
      .def_readwrite("p2", &Phase::p2)
      .def_readwrite("mu", &Phase::mu)
      .def_readwrite("p", &Phase::p)
      .def_readwrite("cdf", &Phase::cdf)
      .def("sample", &Phase::sample);
  py::class_<Material>(m, "Material")
      .def(py::init<>())
      .def_readwrite("mua", &Material::mua)
      .def_readwrite("mus", &Material::mus)
      .def_readwrite("n", &Material::n)
      .def_readwrite("phase", &Material::phase);
  py::class_<VoxelInput>(m, "VoxelInput")
      .def(py::init<>())
      .def_readwrite("nx", &VoxelInput::nx)
      .def_readwrite("ny", &VoxelInput::ny)
      .def_readwrite("nz", &VoxelInput::nz)
      .def_readwrite("dx", &VoxelInput::dx)
      .def_readwrite("dy", &VoxelInput::dy)
      .def_readwrite("dz", &VoxelInput::dz)
      .def_readwrite("labels", &VoxelInput::labels)
      .def_readwrite("materials", &VoxelInput::materials)
      .def_readwrite("n_above", &VoxelInput::n_above)
      .def_readwrite("n_below", &VoxelInput::n_below)
      .def_readwrite("periodic", &VoxelInput::periodic)
      .def_readwrite("cos_incident", &VoxelInput::cos_incident)
      .def_readwrite("seed", &VoxelInput::seed)
      .def_readwrite("packets", &VoxelInput::packets)
      .def_readwrite("w_threshold", &VoxelInput::w_threshold)
      .def_readwrite("roulette_m", &VoxelInput::roulette_m)
      .def_readwrite("record", &VoxelInput::record)
      .def_readwrite("score_voxels", &VoxelInput::score_voxels)
      .def_readwrite("threads", &VoxelInput::threads);
  m.def("uniforms", [](std::uint64_t seed, std::uint64_t id, int n) {
    CounterRandom r(seed, id);
    std::vector<double> v(static_cast<std::size_t>(n));
    for (auto& x : v) x = r.uniform();
    return arr(v);
  });
  auto to_dict = [](const VoxelOutput& o) {
    py::dict d;
    d["sum_r"] = o.sum_r;
    d["sum_r2"] = o.sum_r2;
    d["sum_t"] = o.sum_t;
    d["sum_t2"] = o.sum_t2;
    d["sum_side"] = o.sum_side;
    d["r_specular"] = o.r_specular;
    d["absorbed"] = arr(o.absorbed);
    d["absorbed_voxels"] = arr(o.absorbed_voxels);
    d["id"] = arr(o.id);
    d["side"] = arr(o.side);
    d["x"] = arr(o.x);
    d["y"] = arr(o.y);
    d["ux"] = arr(o.ux);
    d["uy"] = arr(o.uy);
    d["uz"] = arr(o.uz);
    d["w"] = arr(o.w);
    d["zmax"] = arr(o.zmax);
    d["path"] = arr(o.path);
    d["scatterings"] = arr(o.scatterings);
    return d;
  };
  m.def("simulate_voxels", [to_dict](const VoxelInput& in) {
    VoxelOutput o;
    {
      py::gil_scoped_release release;
      o = simulate_voxels(in);
    }
    return to_dict(o);
  });
#ifdef SCATTUCCINO_HAVE_CUDA
  m.attr("HAS_CUDA") = true;
  m.def("cuda_device_count", &cuda_device_count);
  m.def("simulate_voxels_cuda", [to_dict](const VoxelInput& in, int device) {
    VoxelOutput o;
    {
      py::gil_scoped_release release;
      o = simulate_voxels_cuda(in, device);
    }
    return to_dict(o);
  }, py::arg("input"), py::arg("device") = 0);
#else
  m.attr("HAS_CUDA") = false;
  m.def("cuda_device_count", []() { return 0; });
#endif
}
