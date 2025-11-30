#include "../pybind_core.hpp"

template <typename T>
void bind_gemv(py::array_t<T, py::array::forcecast> a,
               py::array_t<T, py::array::forcecast> b,
               py::array_t<T, py::array::forcecast> c);

template <typename T>
py::array_t<T> bind_gemm(py::array_t<T, py::array::forcecast> a,
                         py::array_t<T, py::array::forcecast> b);

template <typename T>
void bind_gemv(py::array_t<T, py::array::forcecast> a,
               py::array_t<T, py::array::forcecast> b,
               py::array_t<T, py::array::forcecast> c);

template <typename T>
py::array_t<T> bind_gemv(py::array_t<T, py::array::forcecast> a,
                         py::array_t<T, py::array::forcecast> b);
