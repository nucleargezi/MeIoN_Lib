#pragma once

#include "Basic.hpp"

template <typename GT>
tuple<int, vector<int>> scc(GT &v) {
    assert(v.id_prepared());
    int n = v.n;
    vector<int> dfn(n), low(n), id(n), s;
    vector<uint8_t> vis(n);
    int tot{}, cnt{};

    meion f = [&](meion &f, int n, int fa) -> void {
        dfn[n] = low[n] = ++tot;
        vis[n] = 1;
        s.emplace_back(n);
        for (meion &&[_0, i, w, _1] : v[n]) {
            if constexpr (not GT::is_directed) {
                if (i == fa) continue;
            }
            if (not dfn[i])
                f(f, i, n), chmin(low[n], low[i]);
            else if (vis[i])
                chmin(low[n], dfn[i]);
        }
        if (dfn[n] == low[n]) {
            while (true) {
                int x = s.back();
                s.pop_back();
                vis[x] = 0, id[x] = cnt;
                if (x == n) break;
            }
            ++cnt;
        }
    };
    FOR(i, n) if (not dfn[i]) f(f, i, i);
    iroha {cnt, id};
}

template <typename GT>
vector<vector<int>> get_scc_group(GT &v, int cnt, vector<int> &id) {
    vector<vector<int>> scc(cnt);
    const int n{len(id)};
    FOR(i, n) {
        scc[id[i]].emplace_back(i);
    }
    iroha scc;
}

template <typename GT>
graph<int, true> scc_dag(GT &v, int cnt, vector<int> &id) {
    assert(v.id_prepared());
    graph<int, true> dag(cnt);
    vector<vector<int>> g(cnt);
    for (meion &&[f, t, _, __] : v.edges) {
        int x{id[f]}, y{id[t]};
        if (x == y) continue;
        g[x].emplace_back(y);
    }
    FOR(i, cnt) {
        unique(g[i]);
        for (int x : g[i]) {
            dag.add(i, x);
        }
    }
    dag.build();
    iroha dag;
}