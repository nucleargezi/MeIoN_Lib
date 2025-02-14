def merge_files(file1_path, file2_path, output_path):
    """
    将两个文件的内容合并到新文件中
    
    参数:
    file1_path (str): 第一个源文件路径
    file2_path (str): 第二个源文件路径
    output_path (str): 合并后的输出文件路径
    """
    try:
        with open(output_path, 'w', encoding='utf-8') as output_file:
            # 写入第一个文件内容
            output_file.write("#### Categorization\n")
            with open(file1_path, 'r', encoding='utf-8') as file:
                output_file.write(file.read())
                output_file.write('\n')  # 添加换行分隔
            output_file.write('\n')
            # 写入第二个文件内容
            with open(file2_path, 'r', encoding='utf-8') as file:
                output_file.write(file.read())
                
        print(f"成功合并文件到 {output_path}")
        
    except FileNotFoundError as e:
        print(f"错误: 文件未找到 - {e.filename}")
    except Exception as e:
        print(f"发生错误: {str(e)}")

if __name__ == "__main__":
    # 在此处修改需要合并的文件名
    FILE1 = "tps.md"    # 第一个源文件名
    FILE2 = "test.md"    # 第二个源文件名
    OUTPUT = "README.md"  # 合并后的文件名
    merge_files(FILE1, FILE2, OUTPUT)