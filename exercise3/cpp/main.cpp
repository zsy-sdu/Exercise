#include "matrix.hpp"
#include <iostream>

int main() {

  auto data = new double[4]{1.0, 1.0, 1.0, 1.0};
  Matrix<double> H(2, 2, data);
  auto [e, c] = eigh(H);
  return 0;
}