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
