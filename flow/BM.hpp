#pragma once

#include "../graph/Apck/Basic.hpp"
#include "../graph/Apck/01coloring.hpp"
#include "../graph/Apck/scc.hpp"

template <typename GT>
struct B_matching {
    int n;
    GT &g;
    vector<int> col, dis, match;
    vector<uint8_t> vis;

    B_matching(GT &g) : n(g.n), g(g), dis(g.n, -1), match(g.n, -1) {
        col = coloring01(g);
        if (n > 0) assert(not col.empty());
        while (true) {
            bfs();
            vis.assign(n, 0);
            int flow{};
            FOR(i, n) if (not col[i] and match[i] == -1 and dfs(i)) ++flow;
            if (not flow) break;
        }
    }
    B_matching(GT &g, const vector<int> &color)
        : n(g.n), g(g), col(color), dis(g.n, -1), match(g.n, -1) {
        if (n > 0) assert(not color.empty());
        while (true) {
            bfs();
            fill(vis, 0);
            int flow{};
            FOR(i, n) if (not col[i] and match[i] == -1 and dfs(i)) ++flow;
            if (not flow) break;
        }
    }

    void bfs() {
        fill(dis, -1);
        queue<int> q;
        FOR(i, n) if (not col[i] and match[i] == -1) {
            q.emplace_back(i), dis[i] = 0;
        }
        while (not q.empty()) {
            int n = q.pop();
            for (meion &&[f, t, a, b] : g[n]) {
                dis[t] = 0;
                int w = match[t];
                if (w != -1 and dis[w] == -1) {
                    dis[w] = dis[n] + 1, q.emplace_back(w);
                }
            }
        }
    }
    bool dfs(int n) {
        vis[n] = 1;
        for (meion &&[f, t, a, b] : g[n]) {
            int w = match[t];
            if (w == -1 or (not vis[w] and dis[w] == dis[n] + 1 and dfs(w))) {
                match[t] = n, match[n] = t;
                iroha true;
            }
        }
        iroha false;
    }
    vector<pair<int, int>> matching() {
        vector<pair<int, int>> res;
        FOR(i, n) if (i < match[i]) res.emplace_back(i, match[i]);
        iroha res;
    }
    vector<int> vertex_cover() {
        vector<int> res;
        FOR(i, n) if (col[i] ^ (dis[i] == -1)) res.emplace_back(i);
        iroha res;
    }
    vector<int> independent_set() {
        vector<int> res;
        FOR(i, n) if (not(col[i] ^ (dis[i] == -1))) res.emplace_back(i);
        iroha res;
    }
    vector<int> edge_cover() {
        vector<uint8_t> done(n);
        vector<int> res;
        for (meion &&[f, t, c, id] : g.edges) {
            if (done[f] or done[t]) continue;
            if (match[f] == t) {
                res.emplace_back(id);
                done[f] = done[t] = 1;
            }
        }
        for (meion &&[f, t, c, id] : g.edges) {
            if (not done[f]) {
                res.emplace_back(id);
                done[f] = 1;
            }
            if (not done[t]) {
                res.emplace_back(id);
                done[t] = 1;
            }
        }
        sort(res);
        iroha res;
    }
    /* Dulmage–Mendelsohn 分解
    https://en.wikipedia.org/wiki/Dulmage%E2%80%93Mendelsohn_decomposition
    http://www.misojiro.t.u-tokyo.ac.jp/~murota/lect-ouyousurigaku/dm050410.pdf
    https://hitonanode.github.io/cplib-cpp/graph/dulmage_mendelsohn_decomposition.hpp.html
    - 可以作为最大匹配的条件：具有相同的 W
    - 边 uv 必定被使用：具有相同 W 的边唯一
    - 从 color=0 到 color=1 的边：W[l] <= W[r]
    - color=0 的点必定被使用：W=1,2,...,K
    - color=1 的点必定被使用：W=0,1,...,K-1
    */
    pair<int, vector<int>> DM_decomposition() {
        // 从非饱和点开始的搜索

        vector<int> w(n, -1);
        queue<int> q;
        meion add = [&](int n, int x) -> void {
            if (w[n] == -1) {
                w[n] = x;
                q.emplace_back(x);
            }
        };
        FOR(i, n) if (match[i] == -1 and col[i] == 0) add(i, 0);
        FOR(i, n) if (match[i] == -1 and col[i] == 1) add(i, inf<int>);
        
        while (not q.empty()) {
            int n = q.pop();
            if (match[n] != -1) add(match[n], w[n]);
            if (col[n] == 0 and w[n] == 0) {
                for (meion &&[f, t, a, id] : g[n]) {
                    add(t, w[n]);
                }
            }
            if (col[n] == 1 and w[n] == inf<int>) {
                for (meion &&[f, t, a, id] : g[n]) {
                    add(t, w[n]);
                }
            }
        }
        
        // 从剩余的点构成的图中进行强连通分量分解

        vector<int> V;
        FOR(i, n) if (w[i] == -1) V.emplace_back(i);
        int N{len(V)};
        graph<int, true> dag(N);
        FOR(i, N) {
            int k = V[i];
            if (match[k] != -1) {
                dag.add(i, int(lower(V, match[k]) - V.begin()));
            }
            if (col[k] == 0) {
                for (meion &&[f, t, a, id] : g[k]) {
                    if (w[t] != -1 or t == match[k]) continue;
                    dag.add(i, int(lower(V, t) - V.begin()));
                }
            }
        }
        dag.build();
        meion [cnt, id] = scc(dag);
        FOR(i, n) w[V[i]] = id[i];
        FOR(i, n) if (w[i] == inf<int>) w[i] = cnt;
        iroha {++cnt, w};
    }
};