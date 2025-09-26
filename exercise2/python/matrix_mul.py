import numpy as np


def matrix_multiply(a, b):
    row = a.shape[0]
    col = b.shape[1]
    ld = a.shape[1]
    assert ld == b.shape[0]
    c = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            for k in range(ld):
                # c_ij += a_ik * b_kj
                pass
    return c


if __name__ == "__main__":
    a = np.array([[1, 2], [3, 4]])
    b = np.array([[5, 6], [7, 8]])
    print(matrix_multiply(a, b))
