import os
import re

def count_cpp_files(directory):
    cpp_count = 0
    for root, _, files in os.walk(directory):
        cpp_count += sum(1 for file in files if file.endswith('.hpp'))
    return cpp_count

def count_cpp_lines_filtered(directory):
    total_lines = 0
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.hpp'):
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        in_main_function = False
                        for line in f:
                            stripped_line = line.strip()
                            # 跳过空行
                            if not stripped_line:
                                continue
                            # 跳过注释行
                            if stripped_line.startswith("//") or stripped_line.startswith("/*") or stripped_line.endswith("*/"):
                                continue
                            # 跳过头文件
                            if stripped_line.startswith("#include"):
                                continue
                            # 跳过宏定义
                            if stripped_line.startswith("#define"):
                                continue
                            # 跳过 main 函数内容
                            if re.match(r"int\s+main\s*\(", stripped_line):  # 匹配 "int main("
                                in_main_function = True
                                continue
                            if in_main_function:
                                if "{" in stripped_line:  # 进入函数体
                                    continue
                                if "}" in stripped_line:  # 离开函数体
                                    in_main_function = False
                                    continue
                                continue
                            
                            # 如果不是以上情况，则计入代码行数
                            total_lines += 1
                except Exception as e:
                    print(f"无法读取文件 {file_path}: {e}")
    return total_lines
if __name__ == "__main__":
    current_directory = os.getcwd()  # 获取当前文件夹路径
    cpp_file_count = count_cpp_files(current_directory)
    print(f"cnt = {cpp_file_count}")

if __name__ == "__main__":
    current_directory = os.getcwd()  # 获取当前文件夹路径
    total_cpp_lines = count_cpp_lines_filtered(current_directory)
    print(f"lines = {total_cpp_lines}")
