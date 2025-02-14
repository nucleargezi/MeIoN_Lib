import os

# 获取当前目录下的所有.hpp文件
def get_files_in_directory(directory):
    files = []
    for root, dirs, filenames in os.walk(directory):
        for filename in filenames:
            # 只选择以.hpp结尾的文件
            if filename.endswith('.hpp'):
                # 获取相对路径
                relative_path = os.path.relpath(os.path.join(root, filename), directory)
                files.append(relative_path)
    return sorted(files)

def write_to_markdown(files, output_filename):
    with open(output_filename, 'w', encoding='utf-8') as md_file:
        # 写入目录部分
        md_file.write('# 目录\n\n')
        for file in files:
            title = file.replace(os.sep, '/')
            md_file.write(f'- [{title}](#{title})\n')
        
        md_file.write('\n')

        # 写入每个文件的内容
        for file in files:
            title = file.replace(os.sep, '/')
            md_file.write(f'## {title}\n\n')
            md_file.write(f'```cpp\n')
            try:
                with open(file, 'r', encoding='utf-8') as f:
                    md_file.write(f.read())
            except Exception as e:
                md_file.write(f"无法读取文件：{file}, 错误信息: {e}")
            md_file.write(f'\n```\n\n')

if __name__ == "__main__":
    # 当前目录
    current_directory = os.getcwd()
    
    # 获取所有.hpp文件
    files = get_files_in_directory(current_directory)
    
    # 将内容写入test.md
    write_to_markdown(files, 'test.md')

    print("生成的markdown文件 'test.md' 已写入当前目录.")
