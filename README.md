![alt text](Z_some_tools/sd_kuk_c01_05.png)
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
#### merge
查询最近一个保存的文件并merge这个文件, 生成一个merge后的文件,  
将文件内容放入剪切板(可以直接粘贴到oj), 删除生成的多余文件  
```bash
python3 MeIoN_Lib/zip_merge.py
```
从学长那里抄的 [Link](https://github.com/YXHXianYu/My-Ancient-Games/tree/main/Shorter%20Machine%20%E8%87%AA%E5%8A%A8%E5%8E%8B%E8%A1%8C%E6%9C%BA)  
查询最近一个保存的文件并merge这个文件, 生成一个merge后的文件,  删除文件内的所有注释,  
将头文件和main函数以外的部分压缩成一行, 将文件内容放入剪切板(可以直接粘贴到oj),  
删除生成的多余文件(用于在有提交行数限制的oj提交)
```bash
python3 MeIoN_Lib/zipzip.py
```
#### to_markdown.py
将文件夹内的所有.hpp后缀的文件内容都按文件夹分块写入一个markdown文档内

#### My_Lib.md / MeIoN's XCPC Library - ICPC2024.pdf
个人XCPC赛时使用的模板, 在比赛结束之后板子还在更新, 仓库里的是最新的

#### flip.sh
对拍用脚本

#### fast.md
指令集加速

#### basic.cpp
初始模板(用于cph插件)