from periodic_table import get_element_info, get_atomic_number, get_atomic_weight, get_element_name

class Molecular:
    def initial_lize(self, file_name):
        print(f"正在从文件 {file_name}初始化分子")

    def __init__(self, file_name):
        self.charge_ = 0   #实例变量，在 __init__ 方法中通过 self. 前缀定义，属于每个对象
        self.spin_ = 1
        self.basis_ = "sto-3g"
        self.atoms_ = []
        self.initial_lize(file_name)#调用实例方法，读取文件

        try:
            with open(file_name, 'r') as f:
                atom_count = int(f.readline().strip())  #每次调用 f.readline() 都会读取文件的下一行
                comment_line = f.readline().strip()
                
                for i in range(atom_count):
                    line = f.readline().strip()
                    if line:
                        parts = line.split() #将字符串按空白字符分割成列表，['O', '0.0', '0.0', '0.0']
                        if len(parts) >= 4:
                            symbol = parts[0]
                            element_info = get_element_info(symbol)
                            
                            atom = {
                                'symbol': symbol,
                                'name': get_element_name(symbol),
                                'atomic_number': get_atomic_number(symbol) or 0,
                                'atomic_weight': get_atomic_weight(symbol) or 0.0,
                                'x': float(parts[1]),
                                'y': float(parts[2]),
                                'z': float(parts[3])
                            }
                            self.atoms_.append(atom)
        except FileNotFoundError:
            print(f"错误: 找不到文件 {file_name}")

    def print(self):
        print("=" * 60)
        print("分子信息")
        print("=" * 60)
        print(f"电荷: {self.charge_}")
        print(f"自旋: {self.spin_}")
        print(f"基组: {self.basis_}")
        print(f"原子数量: {len(self.atoms_)}")
        
        total_weight = sum(atom['atomic_weight'] for atom in self.atoms_)
        print(f"分子量: {total_weight:.3f}")
        
        print("\n原子详细信息:")
        print("-" * 60)
        for i, atom in enumerate(self.atoms_, 1):
            print(f"{i}. {atom['symbol']} ({atom['name']})")
            print(f"   原子序数: {atom['atomic_number']}, 原子量: {atom['atomic_weight']:.3f}")
            print(f"   坐标: ({atom['x']:.6f}, {atom['y']:.6f}, {atom['z']:.6f})")
            print()

    def get_atoms(self):
        return self.atoms_

    def get_charge(self):
        return self.charge_

    def get_spin(self):
        return self.spin_

    def get_basis(self):
        return self.basis_

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