#include "pybind/pybind_core.hpp"

PYBIND11_MODULE(gmev_core, m) {
  m.doc() = "GMEV C++ core GEMM/GEMV bindings";

  // GEMM 自动类型推断
  m.def("gemm", [](py::array a, py::array b) -> py::array {
    if (py::isinstance<py::array_t<float>>(a) &&
        py::isinstance<py::array_t<float>>(b))
      return bind_gemm<float>(a.cast<py::array_t<float>>(),
                              b.cast<py::array_t<float>>());
    else if (py::isinstance<py::array_t<double>>(a) &&
             py::isinstance<py::array_t<double>>(b))
      return bind_gemm<double>(a.cast<py::array_t<double>>(),
                               b.cast<py::array_t<double>>());
    else if (py::isinstance<py::array_t<std::complex<float>>>(a) &&
             py::isinstance<py::array_t<std::complex<float>>>(b))
      return bind_gemm<std::complex<float>>(
          a.cast<py::array_t<std::complex<float>>>(),
          b.cast<py::array_t<std::complex<float>>>());
    else if (py::isinstance<py::array_t<std::complex<double>>>(a) &&
             py::isinstance<py::array_t<std::complex<double>>>(b))
      return bind_gemm<std::complex<double>>(
          a.cast<py::array_t<std::complex<double>>>(),
          b.cast<py::array_t<std::complex<double>>>());
    else
      throw std::runtime_error("gemm: Unsupported dtype");
  });

  // GEMV 自动类型推断
  m.def("gemv", [](py::array a, py::array b) -> py::array {
    if (py::isinstance<py::array_t<float>>(a) &&
        py::isinstance<py::array_t<float>>(b))
      return bind_gemv<float>(a.cast<py::array_t<float>>(),
                              b.cast<py::array_t<float>>());
    else if (py::isinstance<py::array_t<double>>(a) &&
             py::isinstance<py::array_t<double>>(b))
      return bind_gemv<double>(a.cast<py::array_t<double>>(),
                               b.cast<py::array_t<double>>());
    else if (py::isinstance<py::array_t<std::complex<float>>>(a) &&
             py::isinstance<py::array_t<std::complex<float>>>(b))
      return bind_gemv<std::complex<float>>(
          a.cast<py::array_t<std::complex<float>>>(),
          b.cast<py::array_t<std::complex<float>>>());
    else if (py::isinstance<py::array_t<std::complex<double>>>(a) &&
             py::isinstance<py::array_t<std::complex<double>>>(b))
      return bind_gemv<std::complex<double>>(
          a.cast<py::array_t<std::complex<double>>>(),
          b.cast<py::array_t<std::complex<double>>>());
    else
      throw std::runtime_error("gemv: Unsupported dtype");
  });

  py::class_<IntegralWrapper<FunctionType::SPH>>(m, "IntegralSph")
      .def(py::init<py::array_t<int>, py::array_t<int>, py::array_t<double>>())
      .def("intor", &IntegralWrapper<FunctionType::SPH>::intor);

  py::class_<IntegralWrapper<FunctionType::CART>>(m, "IntegralCart")
      .def(py::init<py::array_t<int>, py::array_t<int>, py::array_t<double>>())
      .def("intor", &IntegralWrapper<FunctionType::CART>::intor);

  py::class_<IntegralWrapper<FunctionType::SPINOR>>(m, "IntegralSpinor")
      .def(py::init<py::array_t<int>, py::array_t<int>, py::array_t<double>>())
      .def("intor", &IntegralWrapper<FunctionType::SPINOR>::intor);
}