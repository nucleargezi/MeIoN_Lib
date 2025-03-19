#!/bin/bash

for file in *.out; do
    if [[ $file =~ ^([0-9]+)\.out$ ]]; then
        original_num=${BASH_REMATCH[1]}
        new_num=$((original_num - 1))
        mv -v "$file" "std${new_num}"
    else
        echo "跳过非数字命名的文件: $file"
    fi
done
