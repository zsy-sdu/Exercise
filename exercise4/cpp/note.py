import numpy as np

a = np.array(np.random.rand(2, 3).astype("float32"))
b = np.array(np.random.rand(2, 3).astype("float32"))
c = np.matmul(a, b.T)
print(c)

import gmev_core

d = gmev_core.gemm(a, b.transpose())
print(d)


from pyscf import gto, scf

mol = gto.M(atom="O 0 0 0, H 0 0 1.1, H 1.1 0 0", basis="cc-pvdz")
atm = mol._atm
bas = mol._bas
env = mol._env

driver = gmev_core.IntegralSpinor(
    atm,
    bas,
    env,
)
ovlp = driver.intor("int2e_spinor", "s8")
ovlp_ref = mol.intor("int2e_spinor")

if np.allclose(ovlp, ovlp_ref, rtol=1e-8, atol=1e-10):
    print("ovlp is close")
else:
    print("ovlp is not close")
