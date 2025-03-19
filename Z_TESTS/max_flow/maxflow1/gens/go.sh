#!/bin/bash

for i in {0..14}
do
    if [ $i -lt 5 ]; then
        ./small > "../data/data$i"
    else
        ./large > "../data/data$i"
    fi
    echo "已生成数据文件: data/data$i"
done

echo "全部数据文件生成完成！"
