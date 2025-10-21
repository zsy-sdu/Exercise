import ctypes
from build import build_lib
import numpy as np

# Keep the original _cint and argtypes setup
_cint = build_lib()

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
#restype 是 ctypes 库中的一个属性，用于指定函数的返回值类型。
#这里设置为 ctypes.c_int，表示 _cint1e_ovlp_sph 函数返回一个整数值。

def ovlp_sph(mol):
    atm = mol._atm
    bas = mol._bas
    env = mol._env
    natm = atm.shape[0]
    nbas = bas.shape[0]
    S = np.zeros((nbas, nbas))
    for i in range(nbas):
        for j in range(nbas):
           
            ang_i = bas[i, 1]  # 角动量量子数
            ang_j = bas[j, 1]
            
            di = 2 * ang_i + 1
            dj = 2 * ang_j + 1
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
