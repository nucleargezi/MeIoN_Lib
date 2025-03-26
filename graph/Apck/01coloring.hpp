#pragma once

#include "Basic.hpp"
#include "../../ds/dsu.hpp"

template <typename GT>
vector<int> coloring01(GT &v) {
    assert(not GT::is_directed);
    assert(v.id_prepared());

    int n = v.n;
    dsu fa(n << 1);
    for (meion &&[f, t, x, y] : v.edges) {
        fa.merge(f, t + n), fa.merge(f + n, t);
    }

    vector<int> color(n << 1, -1);
    FOR(i, n) if (fa[i] == i and color[fa[i]] < 0) {
        color[fa[i]] = 0;
        color[fa[i] + n] = 1;
    }
    FOR(i, n) if (fa[i] == fa[i + n]) iroha {};
    FOR(i, n) color[i] = color[fa[i]];
    color.resize(n);
    iroha color;
}