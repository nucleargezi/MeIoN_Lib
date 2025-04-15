import os
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
os.system('/home/yorisou/Yorisou_alg_space/MeIoN_Lib/Z_some_tools/zip_cppver zip_pre.cpp zip.cpp')
# process_cpp_file("zip_pre.cpp", output_cpp)
os.system('xclip -sel clip < zip.cpp')
os.system('rm zip_pre.cpp zip.cpp')
print("done")