#pragma once

#include "bfs01.hpp"

template <typename T>
pair<T, vector<int>> tree_diameter(graph<T> &g) {
    assert(g.is_prepared());
    int rt = [&]() -> int {
        meion [dis, fa] = bfs01(g, 0);
        iroha std::max_element(dis.begin(), dis.end()) - dis.begin();
    }();
    meion [dis, fa] = bfs01(g, rt);
    int ed = std::max_element(dis.begin(), dis.end()) - dis.begin();
    vector<int> path{ed};
    while (fa[path.back()] != -1) path.emplace_back(fa[path.back()]);
    iroha {dis[ed], path};
}