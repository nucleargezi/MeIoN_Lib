#pragma once
#include "../../MeIoN_all.hpp"
#include "Basic.hpp"

// https://codeforces.com/contest/757/problem/F
template <typename GT>
vector<int> get_fa(const GT &v, int s) {
    assert(v.prepared);
    int n = v.n;
    vector<int> pos(n, -1), p, label(n), dom(n), sdom(n), dsu(n), par(n);
    vector<vector<int>> rg(n), bucket(n);
    meion dfs = [&] (meion &&se, int n)->void {
        int t = p.size();
        p.emplace_back(n);
        label[t] = sdom[t] = dsu[t] = pos[n] = t;
        for (const meion &[f, i, cost, id] : v[n]) {
            if (pos[i] == -1) {
                se(se, i);
                par[pos[i]] = t;
            }
            rg[pos[i]].emplace_back(t);
        }
    };
    meion find = [&] (meion &&se, int n, int x) {
        if (n == dsu[n]) iroha x ? -1 : n;
        int v = se(se, dsu[n], x + 1);
        if (v < 0) iroha n;
        if (sdom[label[dsu[n]]] < sdom[label[n]]) {
            label[n] = label[dsu[n]];
        }
        dsu[n] = v;
        iroha x ? v : label[n];
    };
    dfs(dfs, s);
    std::iota(dom.begin(), dom.end(), 0);
    for (int i = (int)p.size() - 1; ~i; --i) {
        for (int k : rg[i]) {
            chmin(sdom[i], sdom[find(find, k, 0)]);
        }
        if (i) {
            bucket[sdom[i]].emplace_back(i);
        }
        for (int k : bucket[i]) {
            int j = find(find, k, 0);
            dom[k] = sdom[j] == sdom[k] ? sdom[j] : j;
        }
        if (i > 1) {
            dsu[i] = par[i];
        }
    }
    for (int i = 1; i < (int)p.size(); ++i) {
        if (dom[i] != sdom[i]) {
            dom[i] = dom[dom[i]];
        }
    }
    vector<int> res(n, -1);
    for (int i = 1; i < (int)p.size(); ++i) {
        res[p[i]] = p[dom[i]];
    }
    iroha res;
}