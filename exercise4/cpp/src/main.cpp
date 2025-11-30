#include "cblas.h"
#include "core/integral.hpp"
#include "core/matmul.hpp"
#include <complex>
#include <iostream>

int main() {
  std::vector<std::complex<float>> a = {
      std::complex<float>(1.0f, 1.0f), std::complex<float>(2.0f, 2.0f),
      std::complex<float>(3.0f, 3.0f), std::complex<float>(4.0f, 4.0f),
      std::complex<float>(5.0f, 5.0f), std::complex<float>(6.0f, 6.0f)};
  std::vector<std::complex<float>> b = {
      std::complex<float>(1.0f, 1.0f), std::complex<float>(2.0f, 2.0f),
      std::complex<float>(3.0f, 3.0f), std::complex<float>(4.0f, 4.0f),
      std::complex<float>(5.0f, 5.0f), std::complex<float>(6.0f, 6.0f)};
  std::vector<std::complex<float>> c(4, std::complex<float>(0.0f, 0.0f));

  core::gemm<std::complex<float>>(CblasColMajor, CblasNoTrans, CblasNoTrans,
                                  a.data(), b.data(), c.data(), 2, 2, 3);

  for (int i = 0; i < c.size(); i++) {
    std::cout << c[i] << " ";
  }
  std::cout << std::endl;

  return 0;
}