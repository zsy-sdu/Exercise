from math import comb
from davidson import davidson_solver


def get_det_num(ncas, neleca, nelecb):
    return comb(ncas, neleca) * comb(ncas, nelecb)


if __name__ == "__main__":
    from pyscf import gto, scf, fci, ao2mo
    import numpy as np

    mol = gto.M(
        atom="H 0 0 0; H 0 0 1",
        basis="6-31g",
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.kernel()

    norb = mol.nao_nr()
    nelec = mol.nelec

    h1 = mf.get_hcore()
    h2 = ao2mo.kernel(mol, mf.mo_coeff, aosym="s1")

    H_fci = fci.direct_spin1.pspace(
        h1, h2, norb, nelec, np=get_det_num(norb, nelec[0], nelec[1])
    )[1]

    e = davidson_solver(
        lambda x: np.dot(H_fci, x), H_fci.diagonal(), n_root=1, max_iter=10
    )

    e_ref, vec_ref = fci.direct_spin1.kernel(h1, h2, norb, nelec, nroots=1)

    # assert np.allclose(e[0], e_ref)
    print(f"e = {e[0]:.12f}, e_ref = {e_ref:.12f}")
