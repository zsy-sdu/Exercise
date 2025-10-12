import numpy as np
import time

class Timer:
    def __enter__(self):
        self.start = time.time()
        return self
    
    def __exit__(self, *args):
        self.end = time.time()
        self.interval = self.end - self.start

def matrix_multiply_naive(a, b):
    """
    纯Python实现的矩阵乘法（三重循环）
    """
    row = a.shape[0]
    col = b.shape[1]
    ld = a.shape[1]
    assert ld == b.shape[0], "矩阵维度不匹配，无法相乘"
    
    c = np.zeros((row, col))
    for i in range(row):
        for j in range(col):
            for k in range(ld):
                c[i, j] += a[i, k] * b[k, j]
    return c

def matrix_multiply_numpy(a, b):
    """
    使用NumPy的dot函数
    """
    return np.dot(a, b)

def matrix_multiply_blas(a, b):
    """
    调用BLAS库进行矩阵乘法
    在NumPy中，dot函数默认会调用BLAS
    """
    return np.dot(a, b)

def compare_performance():
    """
    比较不同方法的性能
    """
    # 创建测试矩阵
    sizes = [10, 50, 100, 200]  # 测试不同大小的矩阵
    
    for size in sizes:
        print(f"\n矩阵大小: {size}x{size}")
        
        # 生成随机矩阵
        a = np.random.rand(size, size)
        b = np.random.rand(size, size)
        
        # 测试纯Python实现
        if size <= 100:  # 对于大矩阵，纯Python太慢，跳过
            with Timer() as t:
                c1 = matrix_multiply_naive(a, b)
            print(f"纯Python实现: {t.interval:.6f} 秒")
        
        # 测试NumPy实现
        with Timer() as t:
            c2 = matrix_multiply_numpy(a, b)
        print(f"NumPy dot函数: {t.interval:.6f} 秒")
        
        # 验证结果一致性（对于小矩阵）
        if size <= 100:
            # 检查结果是否一致（允许浮点数误差）
            diff = np.max(np.abs(c1 - c2))
            print(f"结果差异: {diff:.10f}")

def test_small_matrix():
    """
    测试小矩阵，验证正确性
    """
    print("=== 小矩阵测试 ===")
    a = np.array([[1, 2], [3, 4]])
    b = np.array([[5, 6], [7, 8]])
    
    print("矩阵 A:")
    print(a)
    print("矩阵 B:")
    print(b)
    
    # 使用纯Python实现
    result_naive = matrix_multiply_naive(a, b)
    print("纯Python实现结果:")
    print(result_naive)
    
    # 使用NumPy实现
    result_numpy = matrix_multiply_numpy(a, b)
    print("NumPy dot函数结果:")
    print(result_numpy)
    
    # 验证正确结果应该是 [[19, 22], [43, 50]]
    expected = np.array([[19, 22], [43, 50]])
    print("预期结果:")
    print(expected)

if __name__ == "__main__":
    # 测试小矩阵的正确性
    test_small_matrix()
    
    # 比较性能
    print("\n=== 性能比较 ===")
    compare_performance()
    
    # 额外测试：展示BLAS的优势
    print("\n=== BLAS性能展示 ===")
    large_size = 500
    a_large = np.random.rand(large_size, large_size)
    b_large = np.random.rand(large_size, large_size)
    
    with Timer() as t:
        result_large = matrix_multiply_blas(a_large, b_large)
    print(f"大矩阵 {large_size}x{large_size} 乘法时间: {t.interval:.3f} 秒")