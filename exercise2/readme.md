## 模板元编程初试（C++ Only）

### 编译期运算

在编程中，很多的操作是不需要在程序运行的时候执行的，而是在编译期就可以完成的。例如相对论量化中经常需要计算的一个系数，和光速有关：

```c++
#include <iostream>

int main() {
	// \alpha=1/c^2
  double c = 137.03599967994;
  auto alpha = 1 / (c * c);
  std::cout << alpha << std::endl;
  return 0;
}
```

现在思考，计算alpha的操作一定需要CPU在运行的时候执行指令吗？我们自己摁计算器也可以直接写出alpha的值，所以这个操作显然是没有必要的：

```c++
#include <iostream>

int main() {

  // \alpha=1/c^2
  //   double c = 137.03599967994;
  //   auto alpha = 1 / (c * c);
  auto alpha = 5.32514e-05;
  std::cout << alpha << std::endl;
  return 0;
}
```

但是有时候，我们为了保障可读性，需要 `auto alpha = 1 / (c * c);` ，那么如何实现这个操作呢？只需要给`c`和`alpha`，标记为`constexpr`即可：

```c++
#include <iostream>

int main() {

  // \alpha=1/c^2
  constexpr double c = 137.03599967994;
  constexpr auto alpha = 1 / (c * c);
  std::cout << alpha << std::endl;
  return 0;
}
```

此时我们鼠标悬停到`alpha`上，可以看到其值为：`5.32514e-05`

这个不是非常合适的例子告诉我们，能在编译期执行的事情，就不要放在运行期间。

### 自动生成代码

C/C++是一种静态编译型语言，需要我们制定变量的数据类型。例如我们有一个函数`add`

```c++
#include <iostream>

int add(int x, int y) { return x + y; }
int main() {

  std::cout << add(1, 2) << std::endl;
  return 0;
}
```

可以正常输出3。但如果我的两个变量`x`, `y`是`double`类型:

```c++
#include <iostream>

int add(int x, int y) { return x + y; }
int main() {

  double x = 1.1;
  double y = 2.2;

  std::cout << add(x, y) << std::endl;
  return 0;
}
```

此时也会输出3，导致计算结果错误。经典的解决办法是，使用函数重载，在写一个`double`类型的：

```c++
#include <iostream>

int add(int x, int y) { return x + y; }
double add(double x, double y) { return x + y; }
int main() {

  double x = 1.1;
  double y = 2.2;

  std::cout << add(x, y) << std::endl;
  return 0;
}
```

此时会自动调用第二个函数，正常输出3.3。这样显然太麻烦了！能不能让编译器自己生成对应类型的函数？就需要用到模板template的概念：

```C++
#include <iostream>

template <typename T> T add(T x, T y) { return x + y; }
int main() {

  double x = 1.1;
  double y = 2.2;

  std::cout << add(x, y) << std::endl;
  return 0;
}
```

此时编译器帮我们生成了想要的代码。模板大量出现在STL标准库中。

## 任务1: 计时器Timer

基于模板，写一个模板的计时器，并且基于`if constexpr`实现编译期格式化输出

- [ ] 构造时自动开始计时，析构时输出生命周期的时长
- [ ] 基于std::chrono实现
- [ ] 支持不同的时间单位，如h, min, s, ms, us, ns等

```c++
#pragma once
#ifndef TIMER_HPP
#define TIMER_HPP
#include <chrono>
#include <iostream>

// Timer class for measuring elapsed time
template <typename T = std::chrono::milliseconds> struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Duration = T;

  Clock::time_point start_time_wall;

  Timer() {}
  ~Timer() {
    // Print
  }
  void reset() { /* reset start time */ }

  T elapsed() const {
    // Return during time
  }

  friend std::ostream &operator<<(std::ostream &os, const Timer &timer) {
    // TODO:
  }

  template <typename U> static constexpr std::string_view time_unit_name() {
    if constexpr (std::is_same_v<U, std::chrono::seconds>)
      return "s"; // Support More units
    return "?";
  }
}; // Timer
#endif // TIMER_HPP
```

```c++
// main.cpp
#include "timer.hpp"
#include <thread>

int main() {

  Timer<> timer;
  std::this_thread::sleep_for(std::chrono::seconds(1));
  return 0;
}
```

应当输出：Wall: 1000ms (左右)

## 任务2: 矩阵乘法

量子化学的工作流程简单概括为：推公式 -> 把公式写成代码。我们先从简单的矩阵乘法开始：
$$
C_{ij} = \sum_kA_{ik}B_{kj} 
$$

- [ ] C++可以使用`std::vector<T>`按`Row-Major`存储
- [ ] Python可以使用`numpy.ndarray`
- [ ] 使用python库`time`或上文的`Timer`比较时间
- [ ] 与numpy的`numpy.dot`对比速度
- [ ] 调用一个BLAS（这里以OpenBLAS为例）

```C++
template <typename T>
std::vector<T> matrix_multiply(const std::vector<T> &a, const std::vector<T> &b,
                               const int ld) {
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
```

```python
import numpy as np


def matrix_multiply(a, b):
    row = a.shape[0]
    col = b.shape[1]
    ld = a.shape[1]
    assert ld == b.shape[0]
    c = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            for k in range(ld):
                # c_ij += a_ik * b_kj
                pass
    return c


if __name__ == "__main__":
    a = np.array([[1, 2], [3, 4]])
    b = np.array([[5, 6], [7, 8]])
    print(matrix_multiply(a, b))

```
