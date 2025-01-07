import re
import os

# def remove_comments(code):
#     """去除所有C++注释"""
#     # 去除单行注释
#     code = re.sub(r'//.*', '', code)
#     # 去除多行注释
#     code = re.sub(r'/\*.*?\*/', '', code, flags=re.DOTALL)
#     return code

# def process_cpp_file(input_file, output_file):
#     with open(input_file, 'r', encoding='utf-8') as f:
#         code = f.read()
    
#     # 去除注释
#     code = remove_comments(code)

#     # 分割代码行
#     lines = code.splitlines()

#     # 结果列表
#     result = []
#     buffer = []  # 用于保存压缩行的临时缓冲区
#     is_main_function = False

#     for line in lines:
#         # 去除行首尾空白字符
#         stripped_line = line.strip()

#         # 忽略空行
#         if not stripped_line:
#             continue

#         # 保留头文件和宏定义
#         if stripped_line.startswith('#include') or stripped_line.startswith('#define') or stripped_line.startswith('#ifndef') or stripped_line.startswith('#endif') or stripped_line.startswith('#ifdef'):
#             if buffer:
#                 # 将缓冲区中的代码压缩成一行
#                 result.append(' '.join(buffer))
#                 buffer = []
#             result.append(stripped_line)
        
#         # 检查 main 函数
#         elif 'int main' in stripped_line or stripped_line.startswith('int main'):
#             if buffer:
#                 # 将缓冲区中的代码压缩成一行
#                 result.append(' '.join(buffer))
#                 buffer = []
#             result.append(stripped_line)
#             is_main_function = True

#         elif is_main_function:
#             # main 函数部分，直接保留原格式
#             result.append(line)
        
#         else:
#             # 非 main 函数的部分，压缩代码
#             buffer.append(stripped_line)

#     # 压缩缓冲区中剩余的代码
#     if buffer:
#         result.append(' '.join(buffer))
    
#     # 输出到新的 cpp 文件
#     with open(output_file, 'w', encoding='utf-8') as f:
#         for line in result:
#             f.write(line + '\n')


current_directory = os.getcwd()

files = [f for f in os.listdir(current_directory) if os.path.isfile(os.path.join(current_directory, f))]

if not files:
    print("当前目录下没有文件")
else:
    # find
    latest_file = max(files, key=lambda f: os.path.getctime(os.path.join(current_directory, f)))
    print('latest: ', latest_file)

cpp_merge = 'cpp-merge --output zip_pre.cpp ' + latest_file 
os.system(cpp_merge)
output_cpp = 'zip.cpp'
os.system('/home/guiding_star/MeIoN_FILE/MeIoN_alg_space/MeIoN_Lib/Z_some_tools/zip_cppver zip_pre.cpp zip.cpp')
# process_cpp_file("zip_pre.cpp", output_cpp)
os.system('xclip -sel clip < zip.cpp')
os.system('rm zip_pre.cpp zip.cpp')
print("down")