#pragma once
#include "Basic.hpp"
#include "../../flow/maxflow.hpp"
#include "../../ds/dsu.hpp"

template <typename DAG>
vector<int> dag_path_cover(DAG &v) {
    assert(v.prepared);
    static_assert(DAG::is_directed);
    for (meion&& e : v.edges) {
        assert(e.f < e.to);
    }
    int n = v.n;
    int s = n << 1, t = s | 1;
    MaxFlow<int> FL(t + 1, s, t);
    for (int i{}; i < n; ++i) {
        FL.add(s, i << 1 | 1, 1);
        FL.add(i << 1, t, 1);
        FL.add(i << 1, i << 1 | 1, inf<int>);
    }
    for (meion &&e : v.edges) {
        FL.add(e.f << 1 | 1, e.to << 1, inf<int>);
    }
    FL.flow();
    meion path = FL.path_decomposition();
    dsu g(n);
    for (meion &p : path) {
        int x = p[1], y = p[(int)p.size() - 2];
        g.merge(x >> 1, y >> 1);
    }
    vector<int> ans(n, -1);
    int tot{};
    for (int i{}; i < n; ++i) if (g[i] == i) ans[i] = tot++;
    for (int i{}; i < n; ++i) if (g[i] != i) ans[i] = ans[g[i]];
    iroha ans;
}