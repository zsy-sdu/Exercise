import numpy as np

def get_ovlp(mol):
    return mol.intor('int1e_ovlp_sph')

def get_hcore(mol):
    return mol.intor('int1e_kin_sph') + mol.intor('int1e_nuc_sph')

def get_eri(mol):
    return mol.intor('int2e_sph_sph', aosym='s1')

def get_veff(mol, dm=None, I=None):
    if dm is None:
        dm = make_density(mol) # core guess
    if I is None:
        I = get_eri(mol)
        
    # TODO: J = sum_{kl} I_{ijkl} * dm{kl}
    # TODO: K = sum_{ij} I_{ikjl} * dm{kl}
    
    return J ,K
def get_fock(mol, dm=None, I=None):
    
    # TODO: F = H + 2*J - K
    return F
    
def make_density(mol, fock=None, S=None):
    if S is None:
        S = get_ovlp(mol)
    if fock is None:
        fock = get_hcore(mol) # core guess
    # diagonalize Fock matrix
    # TODO: D = C_occ.T * C_occ
    return D
    
def get_energy_elec(mol, dm=None, fock=None):
    if dm is None:
        dm = make_density(mol) # core guess
    if fock is None:
        fock = get_fock(mol, dm)
    # TODO: E_elec = sum_{ij} (H+F)_{ij} * dm{ij}
    return E_elec

def scf(mol, max_iter=100, tol=1e-6):
    # TODO: SCF loop
    pass

if __name__ == '__main__':
    from pyscf import gto
    
    mol = gto.M(
        atom = 'H 0 0 0; H 0 0 1',
        basis = 'sto-3g',
        charge = 0,
        spin = 0,
    )
    
    scf(mol=mol)