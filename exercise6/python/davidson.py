import numpy as np


def gram_schmidt(xs, freeze_cols=0, lindep=1e-13):
    Y = np.array(xs, copy=True)
    rows, cols = Y.shape

    for i in range(freeze_cols, cols):
        col_i = Y[:, i].copy()

        for j in range(i):
            col_j = Y[:, j]

            denom = np.dot(col_j.conj(), col_j)
            if abs(denom) < 1e-15:
                continue

            coeff = np.dot(col_i.conj(), col_j) / denom
            col_i -= coeff * col_j

        nrm = np.linalg.norm(col_i)
        if nrm > lindep:
            col_i /= nrm

        Y[:, i] = col_i

    return Y


def davidson_solver(
    transformer, diagonal, n_root=1, start_dim=1, max_iter=10, residue_tol=1e-6
):
    assert len(diagonal) == diagonal.shape[0], "diagonal must be a 1D array"
    n_dim = diagonal.shape[0]
    assert start_dim >= n_root, "start_dim must be greater than or equal to n_root"

    # initial guess
    search_space = np.identity(n_dim, dtype=diagonal.dtype)[
        :, :start_dim
    ] + 0.01 * np.ones((n_dim, start_dim), dtype=diagonal.dtype)

    Ab_i = np.empty((n_dim, start_dim), dtype=diagonal.dtype)

    # start iteration
    for iter in range(max_iter):
        # project dim
        M = start_dim + iter * n_root

        orthonormal_subspace = None

        if iter == 0:
            orthonormal_subspace = gram_schmidt(search_space, 0)
            for j in range(start_dim):
                Ab_i[:, j] = transformer(orthonormal_subspace[:, j])
        else:
            new_cols = np.empty((n_dim, n_root), dtype=diagonal.dtype)
            orthonormal_subspace = gram_schmidt(search_space, freeze_cols=M - n_root)
            for j in range(n_root):
                new_cols[:, j] = transformer(orthonormal_subspace[:, M - n_root + j])

            Ab_i = np.hstack((Ab_i, new_cols))

        B = np.dot(orthonormal_subspace.T.conj(), Ab_i)
        eigenvals, eigenvecs = np.linalg.eigh(B)

        xi_n = np.zeros((n_root, n_dim), dtype=diagonal.dtype)
        theta_n = np.zeros(n_root, dtype=diagonal.dtype)
        has_converged = np.zeros(n_root, dtype=bool)

        print(f"davidson diagonalization iter {iter}")

        # eig_pairs[0..n_roots-1]
        for n in range(n_root):
            theta_n[n] = eigenvals[n]  # the n-th eigenvalue

            s = eigenvecs[:, n]
            residue_n = np.dot(Ab_i, s)
            residue_norm = np.linalg.norm(
                residue_n - theta_n[n] * np.dot(orthonormal_subspace, s)
            )

            print(
                f"root {n}: theta = {theta_n[n]:.12f}, |residue| = {residue_norm:.12f}"
            )

            if residue_norm < residue_tol:
                has_converged[n] = True

            # update the search space
            for i in range(n_dim):
                xi_n[n, i] = residue_n[i] / (theta_n[n] - diagonal[i])

            xi_norm = np.linalg.norm(xi_n[n])
            xi_n[n] /= xi_norm

        search_space = np.hstack((orthonormal_subspace, xi_n.T))
        if np.all(has_converged):
            print(f"davidson diagonalization converged in {iter + 1} iterations")
            return theta_n

    raise ValueError(f"davidson diagonalization failed after {iter + 1} iterations")


if __name__ == "__main__":
    from pyscf import gto, scf, fci, ao2mo

    mol = gto.M(
        atom="H 0 0 0; H 0 0 1",
        basis="sto-3g",
        verbose=0,
    )
    mf = scf.RHF(mol)
    mf.kernel()

    norb = 2
    nelec = (1, 1)

    h1 = mf.get_hcore()
    h2 = ao2mo.kernel(mol, mf.mo_coeff, aosym="s1")

    H_fci = fci.direct_spin1.pspace(h1, h2, norb, nelec, np=4)[1]

    transformer = lambda x: np.dot(H_fci, x)

    e = davidson_solver(transformer, H_fci.diagonal())
