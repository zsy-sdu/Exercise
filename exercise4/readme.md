# 获取GTO积分

现代量子化学方法通常建立在分子轨道的理论基础上，对波函数做有限基展开。分子轨道是原子轨道的线性组合，原子轨道可以用Slater Type Function（STF）近似描述，被称为Slater Type Orbital（STO）。但STF的积分不好算，可以选择性质近似的Gaussian Type Function（GTF）做线性组合，拟合成为一个STF，这样的做法得到的轨道称为Gaussian Type Orbital。绝大多数孤立体系量子化学软件都是基于GTO，例如：Gaussian, ORCA, PySCF。仅有极少数软件是基于STO的，例如ADF, BDF-S。

例如1s的STO：
$$
\phi_{1s}^{SF}(\zeta,r-R_A)=(\zeta^3/\pi)^{1/2}e^{-\zeta|r-R_A|}
$$
我们用三个GTF拟合一个Contracted Gaussian Functuon（CGF）
$$
\phi_{1s}^{CGF}(\zeta=1.0,r-R_A)=d_{13}\phi_{1s}^{GF}(\alpha_{13})+d_{23}\phi_{1s}^{GF}(\alpha_{23})+d_{33}\phi_{1s}^{GF}(\alpha_{33})
$$
我们做计算之前指定的基组信息（basis）就是公式2的内容

```
BASIS "ao basis" SPHERICAL PRINT
#BASIS SET: (3s) -> [1s]
H    S
      0.3425250914E+01       0.1543289673E+00
      0.6239137298E+00       0.5353281423E+00
      0.1688554040E+00       0.4446345422E+00
#BASIS SET: (6s,3p) -> [2s,1p]
O    S
      0.1307093214E+03       0.1543289673E+00
      0.2380886605E+02       0.5353281423E+00
      0.6443608313E+01       0.4446345422E+00
O    SP
      0.5033151319E+01      -0.9996722919E-01       0.1559162750E+00
      0.1169596125E+01       0.3995128261E+00       0.6076837186E+00
      0.3803889600E+00       0.7001154689E+00       0.3919573931E+00
END
```

## 调用PYSCF查看基组信息

```python
from pyscf import gto

mol = gto.M(
    atom="H 0 0 0;",
    basis="sto-3g",
    charge=0,
    spin=1,
    verbose=5,
)

```

将输出等级verbose设高一点，可以输出更多信息，尽量解释每一行

## GTF积分库（以Libcint为例）

不同原子轨道之间的积分其实就是一系列公式，可以自行推导但意义不大，可以直接调用现有的积分库。例如PySCF的积分库Libcint

- [ ] 我们最基础的Hartree-Fock需要哪些积分？
- [ ] 计算积分前需要哪些已知信息？
- [ ] Libctin的atm，bas，env的内容都是什么意义？

## 基于PySCF直接调用Libcint计算积分

在PySCF中调用Libcint计算积分非常简单：

```python
ovlp = mol.intor("int1e_ovlp")
```

当然，也可以手写Python代码调用C语言库。在Python文件夹下，我写了一个build.py帮助编译并link到libcint动态库，使用时仅需要：

```python
import ctypes
from build import build_lib

# Keep the original _cint and argtypes setup
_cint = build_lib()
```

此时_cint就是编译好的积分库。C是强类型语言，需要提前声明函数的传入类型和返回类型，对于计算球谐重叠积分的函数，libcint中的定义是：

```C
int cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas,
                    int nbas, double *env);
```

所以我们在代码中声明：

```python
_cint.cint1e_ovlp_sph.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.double, ndim=2),
    (ctypes.c_int * 2),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2),
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.double, ndim=1),
]
_cint.cint1e_ovlp_sph.restype = ctypes.c_int
```

我们从PySCF中拿到atm bas env后，遍历所有shell就能得到我们的积分了：

```python
def ovlp_sph(mol):
    atm = mol._atm
    bas = mol._bas
    env = mol._env
    natm = atm.shape[0]
    nbas = bas.shape[0]
    S = np.zeros((nbas, nbas))
    for i in range(nbas):
        for j in range(nbas):
            di = 1  # TODO:
            dj = 1  # TODO:
            buf = np.empty((di, dj), order="F")
            if (
                _cint.cint1e_ovlp_sph(
                    buf,
                    (ctypes.c_int * 2)(i, j),
                    atm,
                    natm,
                    bas,
                    nbas,
                    env,
                )
                == 0
            ):
                raise RuntimeError("cint1e_ovlp_sph failed")
            S[i, j] = buf[0, 0]
    return S

```

对应之前的H原子，得到的结果已经是1:

```python
if __name__ == "__main__":
    from pyscf import gto

    mol = gto.M(
        atom="""
        H 0.0 0.0 0.0
        """,
        basis="sto-3g",
        charge=0,
        spin=1,
        verbose=5,
    )
    S = ovlp_sph(mol)
    print(S)
```

