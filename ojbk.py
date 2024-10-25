import os

newpage = "<div style=\"page-break-after: always;\"></div>"

def generate_markdown():
    # 创建并打开 README.md 文件
    with open("README.md", "w", encoding="utf-8") as md_file:
        # 写入一级标题
        md_file.write("# Template\n\n")
        
        # 添加 CSS 样式，防止在三级标题前换页
        md_file.write("<style>\n")
        md_file.write("h3 { page-break-before: avoid; }\n")
        md_file.write("</style>\n\n")
        
        # 初始化目录内容
        toc = "## 目录\n\n"
        content = ""
        
        # 遍历当前目录下的所有文件夹
        for root, dirs, files in os.walk("."):
            for dir_name in dirs:
                dir_path = os.path.join(root, dir_name)
                hpp_files = [f for f in os.listdir(dir_path) if f.endswith(".hpp")]
                
                if hpp_files:
                    # 写入目录中的二级标题
                    toc += f"- [{dir_name}](#{dir_name.replace(' ', '-').lower()})\n"
                    
                    # 写入二级标题
                    content += f"\n\n## {dir_name}\n\n"
                    
                    # 遍历文件夹中的所有 .hpp 文件
                    for hpp_file in hpp_files:
                        file_path = os.path.join(dir_path, hpp_file)
                        
                        # 写入目录中的三级标题
                        toc += f"  - [{hpp_file}](#{hpp_file.replace(' ', '-').replace('.', '').lower()})\n"
                        
                        # 写入三级标题
                        content += f"### {hpp_file}\n\n"
                        
                        # 读取 .hpp 文件内容并写入到 Markdown 文件中
                        with open(file_path, "r", encoding="utf-8") as hpp_file_content:
                            content += "```hpp\n"
                            content += hpp_file_content.read()
                            content += "\n```\n"
                    
                    # 添加分页符在每个文件夹内容结束后
                    content += "\n" + newpage + "\n"
        
        # 将目录和内容写入到文件中
        md_file.write(toc + "\n\n" + content)

if __name__ == "__main__":
    generate_markdown()