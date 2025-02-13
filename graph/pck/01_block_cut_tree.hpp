#pragma once
#include "00_basic.hpp"
// [n, n + b_block)
template <typename GT>
graph<int, false> block_cut(GT &g) {
    int n = g.n;
    vector<int> low(n), dfn(n), s;
    vector<uint8_t> vis(n);
    s.reserve(n);
    int nxt = n;
    int k = 0;
    vector<pair<int, int>> edges;
    
    meion dfs = [&](meion &dfs, int n, int fa) -> void {
        s.emplace_back(n);
        low[n] = dfn[n] = k++;
        int child = 0;
        for (meion &&e : g[n]) {
            if (e.to == fa) continue;
            if (not vis[e.to]) {
                ++child;
                int slen = (int)s.size();
                dfs(dfs, e.to, n);
                chmin(low[n], low[e.to]);
                if ((fa == -1 and child > 1) or
                    (fa != -1 and low[e.to] >= dfn[n])) {
                    edges.emplace_back(nxt, n);
                    while ((int)s.size() > slen) {
                        edges.emplace_back(nxt, s.back());
                        s.pop_back();
                    }
                    ++nxt;
                }
            } else {
                chmin(low[n], dfn[e.to]);
            }
        }
    };
    for (int i {}; i < n; ++i) {
        if (not vis[i]) {
            dfs(dfs, i, -1);
            for (meion &&x : s) {
                edges.emplace_back(nxt, x);
            }
            ++nxt;
            s.clear();
        }
    }
    graph<int, false> BCT(nxt);
    for (meion &&[x, y] : edges) BCT.add(x, y);
    BCT.build();
    iroha BCT;
}