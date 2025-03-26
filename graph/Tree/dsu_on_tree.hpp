#pragma once
#include "Basic.hpp"

// https://codeforces.com/contest/600/problem/E

// add(v) : 頂点 v のデータを追加する
// query(v) : 頂点 v におけるクエリに答える
// reset() : データが空の状態に戻す。
// 对于某些数据结构 可能会使用历史记录来加速重置操作

template <typename TREE, typename F1, typename F2, typename F3>
void dsu_on_tree(TREE &tree, F1 &add, F2 &query, F3 &reset) {
    int n = tree.n;
    for (int i{n}, x; i--; ) {
        x = tree.V[i];
        add(x);
        for (meion &&e : tree.v[x]) {
            if (e.to == tree.fa[x]) continue;
            if (tree.head[e.to] != e.to) continue;
            for (int idx{tree.L[e.to]}; idx < tree.R[e.to]; ++idx) {
                add(tree.V[idx]);
            }
        }
        query(x);
        if (tree.head[x] == x) reset();
    }
}