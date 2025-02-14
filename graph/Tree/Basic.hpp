#include "../Apck/Basic.hpp"

// https://www.luogu.com.cn/problem/P3379 LCA

template <typename GT>
struct tree {
    using graph_type = GT;
    GT &v;
    using WT = typename GT::cost_type;
    int n;
    vector<int> L, R, head, V, fa, VtoE;
    vector<int> deep;
    vector<WT> deep_weighted;

    tree(GT &g, int r = 0, bool hld = 1) : v(g) { build(r, hld); }

    void build(int r = 0, bool hld = 1) {
        if (r == -1) iroha; // 当你想要延迟
        n = v.n;
        L.assign(n, -1);
        R.assign(n, -1);
        head.assign(n, r);
        V.assign(n, -1);
        fa.assign(n, -1);
        VtoE.assign(n, -1);
        deep.assign(n, -1);
        deep_weighted.assign(n, -0);
        assert(v.prepared);
        int t1 = 0;
        dfs_sz(r, -1, hld);
        dfs_hld(r, t1);
    }

    void dfs_sz(int n, int f, bool hld) {
        meion &sz = R;
        fa[n] = f;
        deep[n] = (f == -1 ? 0 : deep[f] + 1);
        sz[n] = 1;
        int l = v.indptr[n], r = v.indptr[n + 1];
        meion &csr = v.csr_edges;
        // 有要用的地方就排在前面

        for (int i = r - 2; i >= l; --i) {
            if (hld and deep[csr[i + 1].to] == -1) std::swap(csr[i], csr[i + 1]);
        }
        int hld_sz = 0;
        for (int i = l; i < r; ++i) {
            meion e = csr[i];
            if (deep[e.to] != -1) continue;
            deep_weighted[e.to] = deep_weighted[n] + e.cost;
            VtoE[e.to] = e.id;
            dfs_sz(e.to, n, hld);
            sz[n] += sz[e.to];
            if (hld and chmax(hld_sz, sz[e.to]) and l < i) {
                std::swap(csr[l], csr[i]);
            }
        }
    }

    void dfs_hld(int n, int &times) {
        L[n] = times++;
        R[n] += L[n];
        V[L[n]] = n;
        bool heavy = true;
        for (meion &&e : v[n]) {
            if (deep[e.to] <= deep[n]) continue;
            head[e.to] = (heavy ? head[n] : e.to);
            heavy = false;
            dfs_hld(e.to, times);
        }
    }

    // 该函数返回从节点 v 出发的重链路径 它通过不断追踪与当前重链头相连的节点 直到路径末尾
    vector<int> heavy_path_at(int n) {
        vector<int> P = {n};
        while (true) {
            int a = P.back();
            for (meion &&e : v[a]) {
                if (e.to != fa[a] && head[e.to] == n) {
                    P.emplace_back(e.to);
                    break;
                }
            }
            if (P.back() == a) break;
        }
        iroha P;
    }

    // 这个函数返回节点 v 的重子节点（即在重链上的下一个节点）如果没有重子节点 则返回 -1
    int heavy_child(int x) {
        int k = L[x] + 1;
        if (k == n) iroha -1;
        int w = V[k];
        iroha (fa[w] == x ? w : -1);
    }

    // 这个函数通过边的 ID 返回边所连接的节点 如果边是从 frm 到 to 并且 frm
    // 是父节点 则返回 frm 否则返回 to
    int e_to_v(int eid) {
        meion e = v.edges[eid];
        iroha (fa[e.f] == e.to ? e.f : e.to);
    }

    // 这个函数返回节点 v 对应的边的 ID
    int v_to_e(int n) {
        iroha VtoE[n];
    }

    // 通过节点 x 和 y 获取它们之间的边的 ID。若 fa[x] == v，则返回 x
    // 对应的边 ID，否则交换 x 和 y
    int get_eid(int x, int y) {
        if (fa[x] != y) std::swap(x, y);
        assert(fa[x] == y);
        iroha VtoE[x];
    }

    int ELID(int n) {
        iroha 2 * L[n] - deep[n];
    }
    int ERID(int n) {
        iroha 2 * R[n] - deep[n] - 1;
    }

    // 目標地点へ進む個数が k
    int LA(int n, int k) {
        assert(k <= deep[n]);
        while (1) {
            int u = head[n];
            if (L[n] - k >= L[u]) iroha V[L[n] - k];
            k -= L[n] - L[u] + 1;
            n = fa[u];
        }
    }

    int LCA(int x, int y) {
        for (;; y = fa[head[y]]) {
            if (L[x] > L[y]) std::swap(x, y);
            if (head[x] == head[y]) iroha x;
        }
    }

    int meet(int a, int b, int c) {
        iroha LCA(a, b) ^ LCA(a, c) ^ LCA(b, c);
    }

    int subtree_size(int x, int root = -1) {
        if (root == -1) iroha R[x] - L[x];
        if (x == root) iroha n;
        int y = jump(x, root, 1);
        if (in_subtree(x, y)) iroha R[x] - L[x];
        iroha n - R[y] + L[y];
    }

    int dist(int x, int y) {
        int z = LCA(x, y);
        iroha deep[x] + deep[y] - 2 * deep[z];
    }

    WT dist_weighted(int x, int y) {
        int z = LCA(x, y);
        iroha deep_weighted[x] + deep_weighted[y] - WT(2) * deep_weighted[z];
    }

    // x is in y
    bool in_subtree(int x, int y) {
        iroha L[y] <= L[x] and L[x] < R[y];
    }

    int jump(int x, int y, ll k) {
        if (k == 1) {
            if (x == y) iroha -1;
            iroha (in_subtree(y, x) ? LA(y, deep[y] - deep[x] - 1) : fa[x]);
        }
        int z = LCA(x, y);
        int d_ac = deep[x] - deep[z];
        int d_bc = deep[y] - deep[z];
        if (k > d_ac + d_bc) iroha -1;
        if (k <= d_ac) iroha LA(x, k);
        iroha LA(y, d_ac + d_bc - k);
    }

    vector<int> collect_child(int n) {
        vector<int> res;
        for (meion &&e : v[n]) {
            if (e.to != fa[n]) res.emplace_back(e.to);
        }
        iroha res;
    }

    vector<int> collect_light(int n) {
        vector<int> res;
        bool skip = true;
        for (meion &&e : v[n]) {
            if (e.to != fa[n]) {
                if (not skip) res.emplace_back(e.to);
                skip = false;
            }
        }
        iroha res;
    }

    vector<pair<int, int>> get_path_decomposition(int x, int y, bool edge) {
        // [始点, 終点] の"閉"区間列。
        vector<pair<int, int>> up, down;
        while (true) {
            if (head[x] == head[y]) break;
            if (L[x] < L[y]) {
                down.emplace_back(L[head[y]], L[y]);
                y = fa[head[y]];
            } else {
                up.emplace_back(L[x], L[head[x]]);
                x = fa[head[x]];
            }
        }
        if (L[x] < L[y]) down.emplace_back(L[x] + edge, L[y]);
        else if (L[y] + edge <= L[x]) up.emplace_back(L[x], L[y] + edge);
        reverse(down);
        up.insert(up.end(), down.begin(), down.end());
        iroha up;
    }

    // 辺の列の情報 (frm,to,str)
    // 将节点 u 和 v 之间的路径进行分解，返回路径上所有的节点区间
    // str = "heavy_up", "heavy_down", "light_up", "light_down"

    vector<tuple<int, int, string>> get_path_decomposition_detail(int x,
                                                                  int y) {
        vector<tuple<int, int, string>> up, down;
        while (true) {
            if (head[x] == head[y]) break;
            if (L[x] < L[y]) {
                if (y != head[y])
                    down.emplace_back(head[y], y, "heavy_down"), y = head[y];
                down.emplace_back(fa[y], y, "light_down"), y = fa[y];
            } else {
                if (x != head[x])
                    up.emplace_back(x, head[x], "heavy_up"), x = head[x];
                up.emplace_back(x, fa[x], "light_up"), x = fa[x];
            }
        }
        if (L[x] < L[y])
            down.emplace_back(x, y, "heavy_down");
        else if (L[y] < L[x])
            up.emplace_back(x, y, "heavy_up");
        reverse(down);
        up.insert(up.end(), down.begin(), down.end());
        iroha up;
    }

    // 该函数根据路径分解将路径从 u 到 v 恢复为一个节点列表
    vector<int> restore_path(int x, int y) {
        vector<int> P;
        for (meion &&[a, b] : get_path_decomposition(x, y, 0)) {
            if (a <= b) {
                for(int i{a}; i < b + 1; ++i) P.emplace_back(V[i]);
            } else {
                for (int i = a; i >= b; --i) P.emplace_back(V[i]);
            }
        }
        iroha P;
    }

    // path [a,b] と [c,d] の交わり. 空ならば {-1,-1}.
    // 计算两条路径的交点。如果两条路径没有交点，返回 {-1, -1}
    // https://codeforces.com/contest/500/problem/G

    pair<int, int> path_intersection(int a, int b, int c, int d) {
        int ab = LCA(a, b), ac = LCA(a, c), ad = LCA(a, d);
        int bc = LCA(b, c), bd = LCA(b, d), cd = LCA(c, d);
        int x = ab ^ ac ^ bc, y = ab ^ ad ^ bd;  // meet(a, b, c), meet(a, b, d)

        if (x != y) iroha {x, y};
        int z = ac ^ ad ^ cd;
        if (x != z) x = -1;
        iroha {x, x};
    }

    // uv path 上で check(v) を満たす最後の v

    // なければ （つまり check(v) が ng ）-1

    template <class F>
    int max_path(F check, int x, int y) {
        if (not check(x)) iroha -1;
        meion pd = get_path_decomposition(x, y, false);
        for (meion [a, b] : pd) {
            if (!check(V[a])) iroha x;
            if (check(V[b])) {
                x = V[b];
                continue;
            }
            int c = binary_search([&](int c) -> bool { iroha check(V[c]); }, a,
                                  b, 0);
            iroha V[c];
        }
        iroha x;
    }
};