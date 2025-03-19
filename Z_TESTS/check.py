import subprocess
import os
import shutil

def check(path):
    # 编译程序
    print('compiling~')
    compile_proc = subprocess.run(
        ["g++", "-std=c++20", "-O2", path + "test.cpp", "-o", path + "test"],
        capture_output=True
    )
    if compile_proc.returncode != 0:
        print("Compilation failed:")
        print(compile_proc.stderr.decode())
        return

    os.makedirs(path + "out", exist_ok=True)

    for i in range(100):
        input_file = path + f"data/data{i}"
        output_file = path + f"out/out{i}"
        std_file = path + f"std/std{i}"

        if not os.path.exists(input_file):
            break

        # 运行测试程序（带超时控制）
        try:
            with open(input_file, "rb") as fin, open(output_file, "wb") as fout:
                proc = subprocess.run(
                    [path + "test"],
                    stdin=fin,
                    stdout=fout,
                    timeout=2,
                    stderr=subprocess.PIPE
                )
                
                # 检查运行时错误
                if proc.returncode != 0:
                    print(f"task{i} runtime error (exit code {proc.returncode})")
                    if proc.stderr:
                        print("Error message:", proc.stderr.decode().strip())
                    break

        except subprocess.TimeoutExpired:
            print(f"task{i} time limit exceed")
            break

        # 比较输出结果
        diff_proc = subprocess.run(
            ["diff", "-bB", output_file, std_file],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        if diff_proc.returncode != 0:
            print(f"task{i} wrong answer")
            break
        else:
            print(f"task{i} accept")
    shutil.rmtree(path + 'out')
    os.remove(path + "test")
    print("fin")

def checker(foldername):
    current_dir = os.getcwd()
    
    for entry in os.listdir(current_dir):
        entry_path = os.path.join(current_dir, entry)
        
        # 寻找匹配的测试项
        if os.path.isdir(entry_path) and entry == foldername:
            print(f"\n在 '{entry}' 中找到的测试项目：")
            
            # 获取目标目录下的所有条目
            try:
                sub_entries = os.listdir(entry_path)
            except PermissionError:
                print(f"无法访问目录: {entry_path}")
                continue
                
            # 进行测试
            has_subdir = False
            for sub_entry in sub_entries:
                sub_path = os.path.join(entry_path, sub_entry)
                if os.path.isdir(sub_path):
                    print(f"  - {sub_entry}")
                    check(foldername + '/' + sub_entry + '/')
                    has_subdir = True
            
            if not has_subdir:
                print("  该目录中没有子文件夹")    
            return
    print('没有这个板子')

def Library():
    while True:
        print('Y-Judge> ',end='')
        foldername = input()
        if (foldername.lower() == 'exit') :
            return
        checker(foldername)
    

if __name__ == "__main__":
    print('\nYoisou\'s Library checker')
    Library()
