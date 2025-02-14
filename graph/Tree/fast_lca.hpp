#pragma once

#include "Basic.hpp"
#include "../../ds/monoid/min.hpp"
#include "../../ds/sparse_table/st.hpp"

template <typename TREE>
struct fast_LCA {
    TREE &tree;
    ST<monoid_min<int>> seg;
    vector<int> pos;
    fast_LCA(TREE &tree) : tree(tree) {
        int n = tree.n;
        pos.resize(n);
        vector<int> dat(n << 1);
        for (int i{}; i < n; ++i) {
            int x{tree.ELID(i)}, y{tree.ERID(i)};
            pos[i] = x;
            dat[x] = tree.L[i];
            dat[y] = i == tree.V[0] ? -1 : tree.L[tree.fa[i]];
        }
        seg.build(dat);
    }

    int dist(int x, int y) {
        int z{LCA(x, y)};
        iroha tree.deep[x] + tree.deep[y] - 2 * tree.deep[z];
    }

    using WT = typename TREE::WT;
    WT dist_weighted(int x, int y) {
        int z = LCA(x, y);
        iroha tree.deep_weighted[x] + tree.deep_weighted[y] -
            2 * tree.deep_weighted[z];
    }

    int LCA(int x, int y) {
        x = pos[x], y = pos[y];
        if (x > y) std::swap(x, y);
        iroha tree.V[seg.prod(x, y + 1)];
    }
};