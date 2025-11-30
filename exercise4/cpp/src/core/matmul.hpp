#pragma once
#ifndef GMAT_MATMUL_HPP_
#define GMAT_MATMUL_HPP_

#if defined(_USE_OPENBLAS_)
#include <cblas.h>
#elif defined(_USE_MKL_)
#include <mkl.h>
#elif defined(_USE_ACCELERATE_)
#include <Accelerate/Accelerate.h>
#endif

#include <complex>
#include <type_traits>

namespace core {

inline int infer_lda(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE trans,
                     const int m, const int k) {
  if (layout == CblasRowMajor) {
    return (trans == CblasNoTrans) ? k : m;
  } else { // ColMajor
    return (trans == CblasNoTrans) ? m : k;
  }
}

inline int infer_ldb(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE trans,
                     const int k, const int n) {
  if (layout == CblasRowMajor) {
    return (trans == CblasNoTrans) ? n : k;
  } else { // ColMajor
    return (trans == CblasNoTrans) ? k : n;
  }
}

inline int infer_ldc(const CBLAS_ORDER layout, const int m, const int n) {
  return (layout == CblasRowMajor) ? n : m;
}

template <typename T>
void gemm(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE transa,
          const CBLAS_TRANSPOSE transb, const T *a, const T *b, T *c,
          const int m, const int n, const int k) {
  const int lda = infer_lda(layout, transa, m, k);
  const int ldb = infer_ldb(layout, transb, k, n);
  const int ldc = infer_ldc(layout, m, n);

  if constexpr (std::is_same_v<T, float>) {
    cblas_sgemm(layout, transa, transb, m, n, k, 1.0f, a, lda, b, ldb, 0.0f, c,
                ldc);
  } else if constexpr (std::is_same_v<T, double>) {
    cblas_dgemm(layout, transa, transb, m, n, k, 1.0, a, lda, b, ldb, 0.0, c,
                ldc);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    const std::complex<float> alpha(1.0f, 0.0f);
    const std::complex<float> beta(0.0f, 0.0f);
    cblas_cgemm(layout, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta,
                c, ldc);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double> beta(0.0, 0.0);
    cblas_zgemm(layout, transa, transb, m, n, k, &alpha, a, lda, b, ldb, &beta,
                c, ldc);
  } else {
    static_assert(std::false_type::value,
                  "matmul: Unsupported type for matrix multiplication");
  }
}

template <typename T>
const T *gemm(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE transa,
              const CBLAS_TRANSPOSE transb, const T *a, const T *b, const int m,
              const int n, const int k) {
  T *c = new T[m * n];
  gemm<T>(layout, transa, transb, a, b, c, m, n, k);
  return c;
}

inline int infer_lda_gemv(CBLAS_ORDER layout, CBLAS_TRANSPOSE trans, int m,
                          int n) {
  if (layout == CblasRowMajor) {
    return (trans == CblasNoTrans) ? n : m;
  } else { // ColMajor
    return (trans == CblasNoTrans) ? m : n;
  }
}

template <typename T>
void gemv(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE transa, const T *a,
          const T *b, T *c, const int m, const int n) {
  const int lda = infer_lda_gemv(layout, transa, m, n);

  if constexpr (std::is_same_v<T, float>) {
    cblas_sgemv(layout, transa, m, n, 1.0f, a, lda, b, 1, 0.0f, c, 1);
  } else if constexpr (std::is_same_v<T, double>) {
    cblas_dgemv(layout, transa, m, n, 1.0, a, lda, b, 1, 0.0, c, 1);
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    const std::complex<float> alpha(1.0f, 0.0f);
    const std::complex<float> beta(0.0f, 0.0f);
    cblas_cgemv(layout, transa, m, n, &alpha, a, lda, b, 1, &beta, c, 1);
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    const std::complex<double> alpha(1.0, 0.0);
    const std::complex<double> beta(0.0, 0.0);
    cblas_zgemv(layout, transa, m, n, &alpha, a, lda, b, 1, &beta, c, 1);
  } else {
    static_assert(std::false_type::value,
                  "matmul: Unsupported type for matrix multiplication");
  }
}

template <typename T>
const T *gemv(const CBLAS_ORDER layout, const CBLAS_TRANSPOSE transa,
              const T *a, const T *b, const int m, const int n) {
  T *c = new T[m];
  gemv<T>(layout, transa, a, b, c, m, n);
  return c;
}

} // namespace core

#endif // GMAT_MATMUL_HPP_