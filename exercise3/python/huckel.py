import numpy as np

H = np.array([[1.0, 1.0], [1.0, 1.0]])
e, c = np.linalg.eigh(H)

print(e)
print(c)
