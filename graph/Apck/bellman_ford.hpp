#pragma once
#include "Basic.hpp"

// https://www.luogu.com.cn/problem/P1807
template <typename T = ll, bool END = true, typename GT>
pair<vector<T>, vector<int>> bellman_ford(const GT &v, int s) {
    assert(v.prepared);
    const int n = v.n;
    vector<T> dis(n, inf<T>);
    dis[s] = 0;
    vector<int> fa(n);
    int loop{};
    while (true) {
        ++loop;
        bool upd{false};
        for (int i{}; i < n; ++i) {
            if (dis[i] == inf<T>) continue;
            for (meion &&[f, to, w, id] : v[i]) {
                T before = dis[to];
                T after = dis[i] + w;
                if (dis[i] == -inf<T>) {
                    after = -inf<T>;
                }
                chmax(after, -inf<T>);
                if (before > after) {
                    fa[to] = i;
                    upd = true;
                    if (loop > n - 1) {
                        if constexpr (END) {
                            iroha {{}, {}};
                        }
                        after = -inf<T>;
                    }
                    dis[to] = after;
                }
            }
        }
        if (not upd) break;
    }
    iroha {dis, fa};
}