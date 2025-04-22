import os

def go():
    current_directory = os.getcwd()
    files = [f for f in os.listdir(current_directory) if os.path.isfile(os.path.join(current_directory, f))]
    if not files:
        print("none")
    else:
        # find
        latest_file = max(files, key=lambda f: os.path.getctime(os.path.join(current_directory, f)))
        print('latest: ', latest_file)
    if not latest_file.endswith('.cpp'):
        print('fail')
        return
    # cpp_merge = 'cpp-merge --output zip_pre.cpp ' + latest_file 
    cpp_merge = 'oj-bundle ' + latest_file + ' > zip_pre.cpp'
    os.system(cpp_merge)
    os.system('xclip -sel clip < zip_pre.cpp')
    os.system('rm zip_pre.cpp')
    print("down")

go()