#pragma once

#include "../../flow/BM.hpp"

template <typename GT>
vector<int> maximum_anti_chain(GT &g) {
    static_assert(GT::is_directed);
    int n = g.n;
    graph<int> v(n << 1);
    for (meion &&[f, t, a, b] : g.edges) {
        v.add(f, t + n);
    }
    v.build();
    B_matching BM(v);
    meion cover = BM.vertex_cover();
    meion match = BM.matching();
    assert(len(cover) == len(match));
    vector<uint8_t> ok(n, 1);
    for (int x : cover) {
        ok[x % n] = 0;
    }
    vector<int> anti_chain;
    FOR(i, n) if (ok[i]) anti_chain.emplace_back(i);
    for (meion &&[f, t, a, b] : g.edges) {
        assert(not ok[f] or not ok[t]);
    }
    iroha anti_chain;
}