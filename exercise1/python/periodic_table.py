PERIODIC_TABLE = {
    "H": {"name": "氢", "atomic_number": 1, "atomic_weight": 1.008, "group": 1, "period": 1},
    "He": {"name": "氦", "atomic_number": 2, "atomic_weight": 4.003, "group": 18, "period": 1},
    "Li": {"name": "锂", "atomic_number": 3, "atomic_weight": 6.941, "group": 1, "period": 2},
    "Be": {"name": "铍", "atomic_number": 4, "atomic_weight": 9.012, "group": 2, "period": 2},
    "B": {"name": "硼", "atomic_number": 5, "atomic_weight": 10.81, "group": 13, "period": 2},
    "C": {"name": "碳", "atomic_number": 6, "atomic_weight": 12.011, "group": 14, "period": 2},
    "N": {"name": "氮", "atomic_number": 7, "atomic_weight": 14.007, "group": 15, "period": 2},
    "O": {"name": "氧", "atomic_number": 8, "atomic_weight": 15.999, "group": 16, "period": 2},
    "F": {"name": "氟", "atomic_number": 9, "atomic_weight": 18.998, "group": 17, "period": 2},
    "Ne": {"name": "氖", "atomic_number": 10, "atomic_weight": 20.180, "group": 18, "period": 2},
    "Na": {"name": "钠", "atomic_number": 11, "atomic_weight": 22.990, "group": 1, "period": 3},
    "Mg": {"name": "镁", "atomic_number": 12, "atomic_weight": 24.305, "group": 2, "period": 3},
    "Al": {"name": "铝", "atomic_number": 13, "atomic_weight": 26.982, "group": 13, "period": 3},
    "Si": {"name": "硅", "atomic_number": 14, "atomic_weight": 28.086, "group": 14, "period": 3},
    "P": {"name": "磷", "atomic_number": 15, "atomic_weight": 30.974, "group": 15, "period": 3},
    "S": {"name": "硫", "atomic_number": 16, "atomic_weight": 32.065, "group": 16, "period": 3},
    "Cl": {"name": "氯", "atomic_number": 17, "atomic_weight": 35.453, "group": 17, "period": 3},
    "Ar": {"name": "氩", "atomic_number": 18, "atomic_weight": 39.948, "group": 18, "period": 3},
    "K": {"name": "钾", "atomic_number": 19, "atomic_weight": 39.098, "group": 1, "period": 4},
    "Ca": {"name": "钙", "atomic_number": 20, "atomic_weight": 40.078, "group": 2, "period": 4},
    "Fe": {"name": "铁", "atomic_number": 26, "atomic_weight": 55.845, "group": 8, "period": 4},
    "Cu": {"name": "铜", "atomic_number": 29, "atomic_weight": 63.546, "group": 11, "period": 4},
    "Zn": {"name": "锌", "atomic_number": 30, "atomic_weight": 65.38, "group": 12, "period": 4},
    "Ag": {"name": "银", "atomic_number": 47, "atomic_weight": 107.87, "group": 11, "period": 5},
    "Au": {"name": "金", "atomic_number": 79, "atomic_weight": 196.97, "group": 11, "period": 6}
}

def get_element_info(symbol):
    """根据元素符号获取元素信息"""
    return PERIODIC_TABLE.get(symbol)

def get_atomic_number(symbol):
    """获取原子序数"""
    element = get_element_info(symbol)
    return element["atomic_number"] if element else None

def get_atomic_weight(symbol):
    """获取原子量"""
    element = get_element_info(symbol)
    return element["atomic_weight"] if element else None

def get_element_name(symbol):
    """获取元素名称"""
    element = get_element_info(symbol)
    return element["name"] if element else None



