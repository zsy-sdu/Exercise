import numpy as np
def simple_huckel():
  H = np.array([[1.0, 1.0], [1.0, 1.0]])
  e, c = np.linalg.eigh(H)
  print("\n本征值 (分子轨道能级):")
  print(e)
  print("\n本征向量 (分子轨道系数):")
  print(c)
if __name__ == "__main__":
    simple_huckel()