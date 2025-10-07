#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cblas.h>
#include <iostream>
#include <lapacke.h>
#include <list>
#include <omp.h>
#include <stdexcept>

template <typename T> class Matrix {
public:
  Matrix() : rows_(0), cols_(0), data_(nullptr), order_(CblasRowMajor) {}

  Matrix(int rows, int cols, CBLAS_ORDER order = CblasRowMajor)
      : rows_(rows), cols_(cols), data_(new T[rows * cols]), order_(order) {}
  Matrix(int rows, int cols, T *data)
      : rows_(rows), cols_(cols), data_(data),
        order_(CBLAS_ORDER::CblasRowMajor) {
    if (data == nullptr) {
      throw std::invalid_argument("data is nullptr");
    }
  }
  Matrix(int rows, int cols, std::initializer_list<T> list)
      : rows_(rows), cols_(cols), data_(new T[rows * cols]),
        order_(CblasRowMajor) {
    if (list.size() != rows * cols) {
      throw std::invalid_argument("list size not match");
    }
    std::copy(list.begin(), list.end(), data_);
  }

  ~Matrix() {
    if (data_ != nullptr) {
      delete[] data_;
    }
  }

  Matrix(const Matrix &other) {
    rows_ = other.rows_;
    cols_ = other.cols_;
    order_ = other.order_;
    data_ = new T[rows_ * cols_];
    std::memcpy(data_, other.data_, sizeof(T) * rows_ * cols_);
  }
  Matrix(Matrix &&other) {
    // TODO:
  }

  const int rows() const { return rows_; }
  const int cols() const { return cols_; }
  const T *data() const { return data_; }
  T *data() { return data_; }

  // operator()
  T &operator()(int row, int col) {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
      throw std::out_of_range("Matrix subscript out of range");
    }
    if (order_ == CblasRowMajor) {
      // TODO:
    } else {
      // TODO:
    }
  }
  const T &operator()(int row, int col) const {
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
      throw std::out_of_range("Matrix subscript out of range");
    }
    if (order_ == CblasRowMajor) {
      return data_[row * cols_ + col];
    } else {
      // TODO:
    }
  }

  // operator +
  Matrix<T> operator+(const Matrix<T> &other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::invalid_argument("Matrix size not match");
    }
    Matrix<T> result(rows_, cols_);
    // TODO:
    return result;
  }
  // operator -
  Matrix<T> operator-(const Matrix<T> &other) const {
    if (rows_ != other.rows() || cols_ != other.cols()) {
      throw std::invalid_argument("Matrix size not match");
    }
    Matrix<T> result(rows_, cols_);
    // TODO:
    return result;
  }

  // operator *
  Matrix<T> operator*(const Matrix<T> &other) const {
    if (cols_ != other.rows()) {
      throw std::invalid_argument("Matrix size not match");
    }
    Matrix<T> result(rows_, other.cols());
    // TODO:
    return result;
  }
  Matrix<T> operator*(T scalar) const {
    Matrix<T> result(rows_, cols_);
    // TODO:
    return result;
  }
  // operator *=
  Matrix<T> &operator*=(T scalar) {
    // TODO:
    return *this;
  }

  // operator <<
  friend std::ostream &operator<<(std::ostream &os, const Matrix<T> &mat) {
    for (int i = 0; i < mat.rows_; ++i) {
      for (int j = 0; j < mat.cols_; ++j) {
        os << mat(i, j) << " ";
      }
      os << "\n";
    }
    return os;
  }

private:
  int rows_, cols_;
  T *data_;
  CBLAS_ORDER order_;
};

template <typename T> struct MatrixTypeInfo {
  using RealType = T;
};
template <> struct MatrixTypeInfo<std::complex<float>> {
  using RealType = float;
};
template <> struct MatrixTypeInfo<std::complex<double>> {
  using RealType = double;
};

// eigh
template <typename T>
std::pair<std::vector<typename MatrixTypeInfo<T>::RealType>, Matrix<T>>
eigh(const Matrix<T> &H) {
  using RealType = typename MatrixTypeInfo<T>::RealType;
  if (H.rows() != H.cols()) {
    throw std::invalid_argument("Matrix must be square");
  }
  int n = H.rows();
  std::vector<RealType> w(n);
  Matrix<T> eigvec(H);

  int info = 0;
  if constexpr (std::is_same_v<T, float>) {
    info = LAPACKE_ssyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n,
                          w.data());
  } else if constexpr (std::is_same_v<T, double>) {
    info = LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n,
                          w.data());
  } else if constexpr (std::is_same_v<T, std::complex<float>>) {
    info = LAPACKE_cheevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n,
                          w.data());
  } else if constexpr (std::is_same_v<T, std::complex<double>>) {
    info = LAPACKE_zheevd(LAPACK_ROW_MAJOR, 'V', 'U', n, eigvec.data(), n,
                          w.data());
  } else {
    throw std::runtime_error("Unsupported type for eigh");
  }

  if (info != 0) {
    throw std::runtime_error("LAPACK diagonalization failed");
  }
  return {w, eigvec};
}

template <typename T>
std::pair<std::vector<typename MatrixTypeInfo<T>::RealType>, Matrix<T>>
eigh(const Matrix<T> &H, const Matrix<T> &S) {
  // TODO: *sygev
  return {};
}

#endif // MATRIX_HPP