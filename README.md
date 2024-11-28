![alt text](sd_kuk_c01_05.png)
### guiding star's XCPC_Library
#### qq : 60422310
#### Codeforces : MeIoN_is_UMP45, white_unicorn
#### Luogu : MeIoN
#### Atcoder : MeIoN
##### using cpp_merge to bind cpp file (https://github.com/FastAlien/cpp-merge)
#### basic
```bash
g++ -std=c++2a -DMeIoN Ciallo.cpp -o Ciallo
```
```bash
./Ciallo
```
#### my persional usage
    查询最近一个保存的文件并merge这个文件, 生成一个merge后的文件,  
    将文件内容放入剪切板(可以直接粘贴到oj), 删除生成的多余文件
    Query the most recently saved file, merge this file, generate a merged file, 
    copy the content of the file to the clipboard (which can be directly pasted into the online judge),  
    and delete the unnecessary generated file.
```bash
python3 MeIoN_Lib/zip_merge.py
```
    查询最近一个保存的文件并merge这个文件, 生成一个merge后的文件,  删除文件内的所有注释,  
    将头文件和main函数以外的部分压缩成一行, 将文件内容放入剪切板(可以直接粘贴到oj),  
    删除生成的多余文件(用于在有提交行数限制的oj提交)
    Query the most recently saved file, merge this file, generate a merged file,  
    remove all comments from the file, compress the sections excluding the header files and the main function into a single line,  
    copy the content of the file to the clipboard (which can be directly pasted into the online judge),  
    and delete the unnecessary generated file (for submission to online judges with line count limits).
```bash
python3 MeIoN_Lib/zipzip.py
```
#### to_markdown.py
    将文件夹内的所有.hpp后缀的文件内容都按文件夹分块写入一个markdown文档内
    Write the content of all .hpp files within a folder into a markdown document, organized by folder structure.

#### My_Lib.md / MeIoN's XCPC Library - ICPC2024 - Kunming - guiding-star - 博客园.pdf
    个人XCPC赛时使用的模板
    Template used during personal XCPC competition

#### flip.sh
    对拍用脚本
    Script used for pairing

#### fast.md
    指令集加速
    Instruction set acceleration

#### basic.cpp
    初始模板(用于cph插件)
    Initial template (used for cph)