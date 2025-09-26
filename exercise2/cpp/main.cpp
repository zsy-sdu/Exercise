#include "timer.hpp"
#include <cblas.h>
#include <complex>
#include <iomanip>
#include <iostream>
#include <vector>

template <typename T>
std::vector<T> matrix_multiply(const std::vector<T> &a, const std::vector<T> &b,
                               const int ld) {
  Timer<> timer;
  auto row = a.size() / ld;
  auto col = b.size() / ld;
  std::vector<T> c(row * col, T(0));
  for (auto i = 0; i < row; ++i) {
    for (auto j = 0; j < col; ++j) {
      for (auto k = 0; k < ld; ++k) {
        // c_ij += a_ik * b_kj
      }
    }
  }
  return c;
}

template <typename T>
std::vector<T> matrix_multiply_blas(const std::vector<T> &a,
                                    const std::vector<T> &b, const int ld) {
  Timer<> timer;
  auto row = a.size() / ld;
  auto col = b.size() / ld;
  std::vector<T> c(row * col, T(0));
  if constexpr (std::is_same_v<T, float>) {
    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, row, col, ld, 1.0,
                a.data(), ld, b.data(), ld, 0.0, c.data(), ld);
  } else if constexpr (std::is_same_v<T, double>) {
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
  } else {
    throw std::runtime_error("Unsupported type");
  }
}

template <typename T>
void print_matrix(const std::vector<T> &m, const int row, const int col) {
  if (row * col != m.size()) {
    throw std::runtime_error("Invalid matrix size");
  }
  for (auto i = 0; i < row; ++i) {
    for (auto j = 0; j < col; ++j) {
      T val = m[i * col + j];
      if constexpr (std::is_same_v<T, float>) {
        std::cout << std::fixed << std::setprecision(2) << val << " ";
      }
      // TODO: Add other types
    }
    std::cout << "\n";
  }
}

int main() { return 0; }