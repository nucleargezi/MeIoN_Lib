import re

def remove_comments(code):
    """去除所有C++注释"""
    # 去除单行注释
    code = re.sub(r'//.*', '', code)
    # 去除多行注释
    code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)
    return code

def process_cpp_file(input_file, output_file):
    with open(input_file, 'r', encoding='utf-8') as f:
        code = f.read()
    
    # 去除注释
    code = remove_comments(code)

    # 分割代码行
    lines = code.splitlines()

    # 结果列表
    result = []
    buffer = []  # 用于保存压缩行的临时缓冲区
    is_main_function = False

    for line in lines:
        # 去除行首尾空白字符
        stripped_line = line.strip()

        # 忽略空行
        if not stripped_line:
            continue

        # 保留头文件和宏定义
        if stripped_line.startswith('#include') or stripped_line.startswith('#define') or stripped_line.startswith('#ifndef') or stripped_line.startswith('#endif') or stripped_line.startswith('#ifdef'):
            if buffer:
                # 将缓冲区中的代码压缩成一行
                result.append(' '.join(buffer))
                buffer = []
            result.append(stripped_line)
        
        # 检查 main 函数
        elif 'int main' in stripped_line or stripped_line.startswith('int main'):
            if buffer:
                # 将缓冲区中的代码压缩成一行
                result.append(' '.join(buffer))
                buffer = []
            result.append(stripped_line)
            is_main_function = True

        elif is_main_function:
            # main 函数部分，直接保留原格式
            result.append(line)
        
        else:
            # 非 main 函数的部分，压缩代码
            buffer.append(stripped_line)

    # 压缩缓冲区中剩余的代码
    if buffer:
        result.append(' '.join(buffer))
    
    # 输出到新的 cpp 文件
    with open(output_file, 'w', encoding='utf-8') as f:
        for line in result:
            f.write(line + '\n')

# 使用示例
input_cpp = 'F_Polynomial_Constructioncmpx.cpp'  # 输入的 C++ 文件路径
output_cpp = 'zip.cpp'  # 输出的 C++ 文件路径
process_cpp_file(input_cpp, output_cpp)
print("down")