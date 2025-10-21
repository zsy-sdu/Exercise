from build import build_lib

try:
    libcint = build_lib()
    print("libcint加载成功！")
    print(f"库对象: {libcint}")
except Exception as e:
    print(f"libcint加载失败: {e}")