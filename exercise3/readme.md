# 矩阵类（C++ Only）

数值计算中经常需要使用矩阵操作，不论Python还是C++都有非常经典常用的库，例如Python的Numpy，C++的Eigen3，Armadillo等。我们这里写一个简单的模板类Matrix，实现简单的运算符重载和矩阵对角化（底层调用BLAS/LAPACK）。

```cpp
#pragma once
#ifndef MATRIX_HPP
#define MATRIX_HPP

#endif // MATRIX_HPP
```

这里这几行的意思是，仅编译一次，防止交叉引用。具体含义查一下就可以

## 矩阵的成员变量和构造函数

我们知道，一个矩阵需要有行数，列数，元素；对于代码来说还需要有存储方式，即一行一行存储的行主序（RowMajor， C-Style）和一列一列存储的列主序（ColMajor，Fortran-Style）。Numpy一般默认是RowMajor，Fortran默认ColMajor。这里我们先不考虑转置与否。我们下文都默认是**<u>RowMajor</u>**。存储方式对CPU的缓存命中非常重要！

```cpp
template <typename T> class Matrix {
private:
  int rows_, cols_;
  T *data_;
  CBLAS_ORDER order_;
  
};
```

这里我们用一个裸指针存储元素。我们发现所有的private变量后都有一个_，这是一个约定俗成的代码规范，这样读代码的人一看就知道这个是私有成员。

我们需要写三种构造函数和一个析构函数。我们先看构造函数：

```cpp
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
    // TODO:
  }
  Matrix(Matrix &&other) {
    // TODO:
  }
```

尝试回答这几个构造函数都是在什么时候用的？

## 运算符重载

我们需要重载：`() + - *`，以满足最基本的矩阵运算。

对于`()`，在RowMajor，是一行一行存储，我们要访问第a行第b列的元素，他在数组`data_`中的位置就是`a * cols_ + b`，反之就是`b * rows_ + a`：

```cpp
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
    // TODO:
  } else {
    // TODO:
  }
}
```

对于矩阵加法，我们需要调用BLAS，[参考资料](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-0/blas-routines.html)，这就用到我们之前的if constexpr语法：

```
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
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::invalid_argument("Matrix size not match");
  }
  Matrix<T> result(rows_, cols_);
  // TODO:
  return result;
}
```

乘法是同理的，我们实现两个，一个是与浮点数乘，一个是矩阵乘法：

```cpp
// operator *
Matrix<T> operator*(const Matrix<T> &other) const {
  if (cols_ != other.rows_) {
    throw std::invalid_argument("Matrix size not match");
  }
  Matrix<T> result(rows_, other.cols_);
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
```

## 矩阵对角化

分为两种，HC=eC和HC=eSC，这个等大家学会Hartree-Fock就知道了。我们两种都写一下，都写在Matrix的类外面。首先在成员外添加类型匹配的模板：

```cpp
template <typename T> struct MatrixTypeInfo_ {
  using RealType = T;
};
template <> struct MatrixTypeInfo_<std::complex<float>> {
  using RealType = float;
};
template <> struct MatrixTypeInfo_<std::complex<double>> {
  using RealType = double;
};
```

这是因为，如果是complex类型的hermit矩阵对角化，其本征值一定是实数，我们的返回值使用`std::pair`封装：

```cpp
// eigh
template <typename T>
std::pair<std::vector<typename MatrixTypeInfo<T>::RealType>, Matrix<T>>
eigh(const Matrix<T> &H) {
  return {};
}

template <typename T>
std::pair<std::vector<typename MatrixTypeInfo<T>::RealType>, Matrix<T>>
eigh(const Matrix<T> &H, const Matrix<T> &S) {
  return {};
}
```

# 第三方矩阵库

- [ ] 学习`Numpy`
- [ ] 熟悉`Numpy`中的n-dimension array的基本算法（矩阵乘法，矩阵对角化）
- [ ] 熟悉`np.einsum`的使用（很重要！）
- [ ] （必须学C++的同学）尝试使用Eigen3（建议），Armadillo或其他应用广泛的矩阵库中的至少一个

# Huckel分子轨道

实现一个Huckel分子轨道计算的程序，可以基于自己上面写好的矩阵库，也可以基于Eigen3等第三方库。

1. 构造Hamiltonian
2. 对角化

```cpp
#include "matrix.hpp"
#include <iostream>

int main() {

  auto data = new double[4]{1.0, 1.0, 1.0, 1.0};
  Matrix<double> H(2, 2, data);
  auto [e, c] = eigh(H);
  return 0;
}
```

```python
import numpy as np

H = np.array([[1.0, 1.0], [1.0, 1.0]])
e, c = np.linalg.eigh(H)

print(e)
print(c)

```

