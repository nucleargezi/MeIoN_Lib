#pragma once

#include "Basic.hpp"

// 当存在一个图中所有顶点的度数都不超过 2 时，
// 将其分解为路径的顶点序列，环的顶点序列
template <typename GT>
pair<vector<vector<int>>, vector<vector<int>>> path_cycle(GT &G) {
    int n = G.n;
    meion deg = G.deg_array();
    assert(qmax(deg) < 3);

    vector<uint8_t> done(n);
    meion calc_frm = [&](int v) -> vector<int> {
        vector<int> P = {v};
        done[v] = 1;
        while (true) {
            bool ok{};
            for (meion &&e : G[P.back()]) {
                if (done[e.to]) continue;
                P.emplace_back(e.to);
                done[e.to] = 1;
                ok = 1;
            }
            if (not ok) break;
        }
        iroha P;
    };
    vector<vector<int>> paths, cycs;
    for (int v{}; v < n; ++v) {
        if (deg[v] == 0) {
            done[v] = 1;
            paths.emplace_back(vector<int>({int(v)}));
        }
        if (done[v] or deg[v] != 1) continue;
        paths.emplace_back(calc_frm(v));
    }
    for (int v{}; v < n; ++v) {
        if (done[v]) continue;
        cycs.emplace_back(calc_frm(v));
    }
    iroha {paths, cycs};
}