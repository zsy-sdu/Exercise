#pragma once
#ifndef GMAT_PYBIND_CORE_HPP_
#define GMAT_PYBIND_CORE_HPP_

#include <omp.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl_bind.h>

#include "../core/integral.hpp"
#include "../core/matmul.hpp"

namespace py = pybind11;
using namespace core;

// Helper: choose layout preference from optional arrays: prefer c (if given),
// then a, then b
inline CBLAS_ORDER choose_layout_opt(const py::array &a,
                                     const py::array *b = nullptr,
                                     const py::array *c = nullptr) {
  auto is_c = [](const py::array &arr) {
    return (arr.flags() & py::array::c_style) != 0;
  };
  auto is_f = [](const py::array &arr) {
    return (arr.flags() & py::array::f_style) != 0;
  };

  if (c) {
    if (is_c(*c))
      return CblasRowMajor;
    if (is_f(*c))
      return CblasColMajor;
  }
  if (is_c(a))
    return CblasRowMajor;
  if (is_f(a))
    return CblasColMajor;
  if (b) {
    if (is_c(*b))
      return CblasRowMajor;
    if (is_f(*b))
      return CblasColMajor;
  }
  // fallback: prefer RowMajor
  return CblasRowMajor;
}

// Helper: check whether array is exactly contiguous in given layout without
// transpose For shape (rows, cols) and given layout, "not transposed" means:
//   RowMajor contiguous: strides == (cols*sizeof, sizeof)
//   ColMajor contiguous: strides == (sizeof, rows*sizeof)
inline bool is_contiguous_nontransposed(const py::array &arr,
                                        CBLAS_ORDER layout) {
  ssize_t rows = arr.shape(0);
  ssize_t cols = arr.shape(1);
  ssize_t s0 = arr.strides(0);
  ssize_t s1 = arr.strides(1);
  ssize_t elem = static_cast<ssize_t>(arr.itemsize());

  if (layout == CblasRowMajor) {
    return (s0 == cols * elem) && (s1 == elem);
  } else { // ColMajor
    return (s0 == elem) && (s1 == rows * elem);
  }
}

// Helper: check whether array is exactly the transpose view w.r.t. chosen
// layout i.e. its memory layout equals the transpose of non-transposed layout
inline bool is_contiguous_transposed_view(const py::array &arr,
                                          CBLAS_ORDER layout) {
  ssize_t rows = arr.shape(0);
  ssize_t cols = arr.shape(1);
  ssize_t s0 = arr.strides(0);
  ssize_t s1 = arr.strides(1);
  ssize_t elem = static_cast<ssize_t>(arr.itemsize());

  if (layout == CblasRowMajor) {
    // non-transposed expected (rows,cols) -> (cols*elem, elem)
    // transposed view would have (elem, rows*elem)
    return (s0 == elem) && (s1 == rows * elem);
  } else { // ColMajor
    // non-transposed expected (elem, rows*elem)
    // transposed view would have (cols*elem, elem)
    return (s0 == cols * elem) && (s1 == elem);
  }
}

// Determine trans flag for a given array (shape known) relative to chosen
// layout. If array matches non-transposed stride pattern => NoTrans. If array
// matches transposed-view stride pattern => Trans. Otherwise, try to detect
// "general" case: if innermost stride equals sizeof(T) then treat as NoTrans,
// else if outermost stride equals sizeof(T) treat as Trans. (last-resort
// heuristic)
inline CBLAS_TRANSPOSE detect_transpose(const py::array &arr,
                                        CBLAS_ORDER layout) {
  if (is_contiguous_nontransposed(arr, layout))
    return CblasNoTrans;
  if (is_contiguous_transposed_view(arr, layout))
    return CblasTrans;

  // Last-resort heuristic: check which axis has stride == itemsize
  ssize_t s0 = arr.strides(0);
  ssize_t s1 = arr.strides(1);
  ssize_t elem = static_cast<ssize_t>(arr.itemsize());
  if (layout == CblasRowMajor) {
    // expected inner stride == elem (axis 1)
    if (s1 == elem && s0 != elem)
      return CblasNoTrans;
    if (s0 == elem && s1 != elem)
      return CblasTrans;
  } else {
    // ColMajor: expected inner stride == elem (axis 0)
    if (s0 == elem && s1 != elem)
      return CblasNoTrans;
    if (s1 == elem && s0 != elem)
      return CblasTrans;
  }
  // fallback
  return CblasNoTrans;
}

////////////////////
// GEMM (in-place c provided)
template <typename T>
void bind_gemm(py::array_t<T, py::array::forcecast> a,
               py::array_t<T, py::array::forcecast> b,
               py::array_t<T, py::array::forcecast> c) {
  if (a.ndim() != 2 || b.ndim() != 2 || c.ndim() != 2)
    throw std::runtime_error("Input arrays must be 2D");

  ssize_t m = a.shape(0);
  ssize_t k1 = a.shape(1);
  ssize_t k2 = b.shape(0);
  ssize_t n = b.shape(1);

  if (k1 != k2)
    throw std::runtime_error("Inner dimensions do not match for GEMM");
  if (c.shape(0) != m || c.shape(1) != n)
    throw std::runtime_error("Output array shape does not match GEMM");

  // Choose layout: prefer c (user-provided), then a, then b
  CBLAS_ORDER layout = choose_layout_opt(a, &b, &c);

  // Deduce transpose flags robustly
  CBLAS_TRANSPOSE transA = detect_transpose(a, layout);
  CBLAS_TRANSPOSE transB = detect_transpose(b, layout);

  auto a_ptr = a.data();
  auto b_ptr = b.data();
  auto c_ptr = c.mutable_data();

  core::gemm<T>(layout, transA, transB, a_ptr, b_ptr, c_ptr,
                static_cast<int>(m), static_cast<int>(n), static_cast<int>(k1));
}

////////////////////
// GEMM (returns new c)
template <typename T>
py::array_t<T> bind_gemm(py::array_t<T, py::array::forcecast> a,
                         py::array_t<T, py::array::forcecast> b) {
  if (a.ndim() != 2 || b.ndim() != 2)
    throw std::runtime_error("Input arrays must be 2D");

  ssize_t m = a.shape(0);
  ssize_t k1 = a.shape(1);
  ssize_t k2 = b.shape(0);
  ssize_t n = b.shape(1);

  if (k1 != k2)
    throw std::runtime_error("Inner dimensions do not match for GEMM");

  // Choose layout: prefer a, then b
  CBLAS_ORDER layout = choose_layout_opt(a, &b, nullptr);

  // Deduce transpose flags
  CBLAS_TRANSPOSE transA = detect_transpose(a, layout);
  CBLAS_TRANSPOSE transB = detect_transpose(b, layout);

  // Allocate output c matching chosen layout:
  py::array_t<T> c;
  if (layout == CblasRowMajor) {
    // C-contiguous: shape (m,n), row major strides (n*sizeof, sizeof)
    c = py::array_t<T>({m, n});
  } else {
    // Fortran-contiguous: provide explicit strides to create F-contiguous numpy
    // array
    ssize_t elem = static_cast<ssize_t>(sizeof(T));
    std::vector<ssize_t> strides = {elem,
                                    elem * m}; // strides in bytes for (m,n)
    py::buffer_info info(nullptr, elem, py::format_descriptor<T>::value, 2,
                         {m, n}, strides);
    // create empty array with given shape & strides by allocating buffer
    // ourselves:
    T *buf = new T[m * n];
    memset(buf, 0, sizeof(T) * m * n);
    auto capsule =
        py::capsule(buf, [](void *p) { delete[] reinterpret_cast<T *>(p); });
    c = py::array_t<T>(py::buffer_info(buf, elem,
                                       py::format_descriptor<T>::value, 2,
                                       {m, n}, strides),
                       capsule);
  }

  auto c_ptr = c.mutable_data();

  core::gemm<T>(layout, transA, transB, a.data(), b.data(), c_ptr,
                static_cast<int>(m), static_cast<int>(n), static_cast<int>(k1));
  return c;
}

////////////////////
// GEMV (in-place c provided)
template <typename T>
void bind_gemv(py::array_t<T, py::array::forcecast> a,
               py::array_t<T, py::array::forcecast> b,
               py::array_t<T, py::array::forcecast> c) {
  if (a.ndim() != 2 || b.ndim() != 1 || c.ndim() != 1)
    throw std::runtime_error("Input arrays must be 2D and 1D");
  ssize_t m = a.shape(0);
  ssize_t n = a.shape(1);
  ssize_t k = b.shape(0);

  if (n != k)
    throw std::runtime_error("Inner dimensions do not match for GEMV");
  if (c.shape(0) != m)
    throw std::runtime_error("Output array shape does not match GEMV");

  // choose layout prefer c then a
  CBLAS_ORDER layout = choose_layout_opt(a, nullptr, &c);

  CBLAS_TRANSPOSE transA = detect_transpose(a, layout);

  core::gemv<T>(layout, transA, a.data(), b.data(), c.mutable_data(),
                static_cast<int>(m), static_cast<int>(n));
}

////////////////////
// GEMV (returns new c)
template <typename T>
py::array_t<T> bind_gemv(py::array_t<T, py::array::forcecast> a,
                         py::array_t<T, py::array::forcecast> b) {
  if (a.ndim() != 2 || b.ndim() != 1)
    throw std::runtime_error("Input arrays must be 2D and 1D");
  ssize_t m = a.shape(0);
  ssize_t n = a.shape(1);
  ssize_t k = b.shape(0);

  if (n != k)
    throw std::runtime_error("Inner dimensions do not match for GEMV");

  CBLAS_ORDER layout = choose_layout_opt(a, &b, nullptr);
  CBLAS_TRANSPOSE transA = detect_transpose(a, layout);

  py::array_t<T> result;
  if (layout == CblasRowMajor) {
    result = py::array_t<T>({m});
  } else {
    // vector layout is same in memory, just allocate C-contiguous vector
    result = py::array_t<T>({m});
  }

  core::gemv<T>(layout, transA, a.data(), b.data(), result.mutable_data(),
                static_cast<int>(m), static_cast<int>(n));
  return result;
}

template <FunctionType FT>
  requires(FT == FunctionType::CART || FT == FunctionType::SPH ||
           FT == FunctionType::SPINOR)
class IntegralWrapper {
public:
  Integral<FT> impl;
  using T = typename Integral<FT>::T;

  // Constructor
  IntegralWrapper(py::array_t<int> atm, py::array_t<int> bas,
                  py::array_t<double> env)
      : impl(atm.shape(0), atm.mutable_data(), bas.shape(0), bas.mutable_data(),
             env.mutable_data(), static_cast<int>(env.size())) {}

  py::array_t<T> intor(const std::string &name, const std::string &sym = "s1") {
    if (name.find("1e") != std::string::npos) {
      if (name == ovlp_key()) {
        return wrap_int1e(impl.ovlp());
      } else if (name == kin_key()) {
        return wrap_int1e(impl.kin());
      } else if (name == nuc_key()) {
        return wrap_int1e(impl.nuc());
      } else
        throw std::runtime_error("Unknown integral name: " + name);

    } else if (name.find("2e") != std::string::npos) {
      if (sym == "s8" && name == twoe_key()) {
        return wrap_int2e(impl.int2e_s8());
      } else if (sym == "s4" && name == twoe_key()) {
        return wrap_int2e(impl.int2e_s4());
      } else if (sym == "s1" && name == twoe_key()) {
        return wrap_int2e(impl.int2e_s1());
      } else
        throw std::runtime_error("Unknown integral name: " + name);
    }
    throw std::runtime_error("Unknown integral name: " + name);
  }

private:
  static inline std::string ovlp_key() {
    if constexpr (FT == FunctionType::SPH)
      return "int1e_ovlp_sph";
    if constexpr (FT == FunctionType::CART)
      return "int1e_ovlp_cart";
    if constexpr (FT == FunctionType::SPINOR)
      return "int1e_ovlp_spinor";
  }
  static inline std::string kin_key() {
    if constexpr (FT == FunctionType::SPH)
      return "int1e_kin_sph";
    if constexpr (FT == FunctionType::CART)
      return "int1e_kin_cart";
    if constexpr (FT == FunctionType::SPINOR)
      return "int1e_kin_spinor";
  }
  static inline std::string nuc_key() {
    if constexpr (FT == FunctionType::SPH)
      return "int1e_nuc_sph";
    if constexpr (FT == FunctionType::CART)
      return "int1e_nuc_cart";
    if constexpr (FT == FunctionType::SPINOR)
      return "int1e_nuc_spinor";
  }
  static inline std::string twoe_key() {
    if constexpr (FT == FunctionType::SPH)
      return "int2e_sph";
    if constexpr (FT == FunctionType::CART)
      return "int2e_cart";
    if constexpr (FT == FunctionType::SPINOR)
      return "int2e_spinor";
  }

  py::array_t<T> wrap_int1e(T *buf) {
    ssize_t nao = impl.nao_;
    auto capsule =
        py::capsule(buf, [](void *p) { delete[] reinterpret_cast<T *>(p); });
    return py::array_t<T>({nao, nao}, {sizeof(T) * nao, sizeof(T)}, buf,
                          capsule);
  }

  py::array_t<T> wrap_int2e(T *buf) {
    ssize_t nao = impl.nao_;
    auto capsule =
        py::capsule(buf, [](void *p) { delete[] reinterpret_cast<T *>(p); });
    return py::array_t<T>({nao, nao, nao, nao},
                          {sizeof(T) * nao * nao * nao, sizeof(T) * nao * nao,
                           sizeof(T) * nao, sizeof(T)},
                          buf, capsule);
  }
};

#endif // GMAT_PYBIND_CORE_HPP_