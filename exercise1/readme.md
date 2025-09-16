# 面向对象编程

程序员，尤其是C++程序员的工作，可以理解为创建一个类，测试这个类并上传。对于电子结构代码，最基本的类就是分子（Molecular），他至少需要包括：原子类型与原子坐标，电荷，自旋多重度（2S+1或2S），基组信息，核-核排斥势。更进一步需要包括，总电子数，总轨道数，基函数的指数/系数项...我们先写一个简单的分子类，暂时不考虑基组。

## 任务

- [ ] 学习C++/Python面向对象编程基础知识
- [ ] 使用一定的数据结构保存元素周期表，方便后续操作
- [ ] 补全空缺的代码

## 注意

- 我们一般输入文件默认的原子坐标单位是Å（1e-10米)，但公式一般的单位是原子单位，注意转换

- xyz格式中没有charge 和 spin的信息，可以定义第二行注释格式，例如

  ```
  3
  0 1 # charge=0 2S=1
  O 0   0     0
  H 0  -0.757 0.587
  H 0   0.757 0.587
  ```

## 私有成员变量

使用一定的数据结构存储原子信息：

```c++
struct Atom {
  std::string type; // element type
  std::size_t Z;    // atomic number
  double x, y, z;   // coordinates
};
```

因为我们有一系列的原子要存，所以可以使用vector：

```c++
private:
  std::vector<Atom> atoms;
```

同时需要分子的电荷量，自旋，基组等信息：

```c++
  std::vector<Atom> atoms_; // atoms in the molecule
  int charge_;              // charge of the molecule
  int spin_;                // spin multiplicity of the molecule
  std::string basis_;       // basis set name
```

## 公共接口

我们需要一个有参的构造函数，传入一个文件名，要求是xyz格式的，然后解析这个文件内容，此时对我们的成员变量做初始化：

```c++
public:
  // Default constructor
  Molecular() = default;
  // Constructor with filename
  explicit Molecular(const std::string &filename);
  // Copy constructor
  Molecular(const Molecular &other);
  // Move constructor
  Molecular(Molecular &&other);
  // Destructor
  ~Molecular() = default;
```

函数的声明在`molecular.h`中，实现在`molecular.cpp`中：

```c++
#include "molecular.h"

Molecular::Molecular(const std::string &filename) {
  /*
  TODO: read filename and parse the file
  */
}
```

此外，请自行学习**拷贝构造函数** **移动构造函数** **析构函数**，并将相关的代码补齐。

我们需要为我们的类写接口，这样才能在外部访问或修改私有成员变量：

```c++
// Getters
const std::vector<Atom> get_atoms() const;
const int get_charge() const;
const int get_spin() const;
const std::string get_basis() const;
// Setters
void set_atoms(const std::vector<Atom> &atoms);
void set_charge(const int charge);
void set_spin(const int spin);
void set_basis(const std::string &basis);
```

同时通过运算符重载，输出成员信息，方便debug：

```c++
// Overload << operator
friend std::ostream &operator<<(std::ostream &os, const Molecular &molecule);
```

## 如何使用？

```c++
#include "molecular.h"
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " filename" << std::endl;
    return 1;
  }
  /*
  const char *filename = argv[1];
  Molecular molecule(filename);
  std::cout << molecule << std::endl;
  */

  return 0;
}
```

## Python版本

```python
class Molecular:
    def __init__(self, file_name):
        # TODO:
        self.charge_ = 0
        self.spin_ = 0
        self.basis_ = None
        self.atoms_ = None
        self.initial_lize()

    def initial_lize(self):
        # TODO:
        pass

    # print the information of the molecule
    def print(self):
        pass

    # Getter
    def get_atoms(self):
        return self.atoms_

    def get_charge(self):
        return self.charge_

    def get_spin(self):
        return self.spin_

    def get_basis(self):
        return self.basis_

    # Setter
    def set_atoms(self, atoms):
        self.atoms_ = atoms

    def set_charge(self, charge):
        self.charge_ = charge

    def set_spin(self, spin):
        self.spin_ = spin

    def set_basis(self, basis):
        self.basis_ = basis


if __name__ == "__main__":
    mol = Molecular("h2o.xyz")
    mol.print()

```

