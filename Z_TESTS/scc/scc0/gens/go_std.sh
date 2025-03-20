#!/bin/bash

for i in {0..14}
do
    ./tt < "data/data$i" > "std/std$i"
    echo "已生成数据文件: std/std$i"
done

echo "全部数据文件生成完成！"
