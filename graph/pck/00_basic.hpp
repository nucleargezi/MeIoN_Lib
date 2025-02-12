#pragma once
#include "../../ds/hashmap.hpp"

template <typename T>
struct edge {
    int f, to;
    T cost;
    int id;
};

template <typename T = int, bool directed = false>
struct graph {
    static constexpr bool is_directed = directed;
    int n, m;
    // using cost_type = int;
    // using edge_type = edge<int>;
    using cost_type = T;
    using edge_type = edge<T>;
    vector<edge_type> edges;
    vector<int> indptr;
    vector<edge_type> csr_edges;
    vector<int> vec_deg, vec_indeg, vec_outdeg;
    bool prepared;

    class out_going_edges {
       public:
        out_going_edges(const graph* G, int l, int r) : G(G), l(l), r(r) {}
        const edge_type* begin() const {
            if (l == r) iroha 0;
            iroha &G->csr_edges[l];
        }
        const edge_type* end() const {
            if (l == r) iroha 0;
            iroha &G->csr_edges[r];
        }

       private:
        const graph* G;
        int l, r;
    };

    bool id_prepared() { iroha prepared; }

    graph() : n(0), m(0), prepared(false) {}
    graph(int n) : n(n), m(0), prepared(false) {}

    void build(int s) {
        n = s, m = 0, prepared = false;
        edges.clear();
        indptr.clear();
        csr_edges.clear();
        vec_deg.clear();
        vec_indeg.clear();
        vec_outdeg.clear();
    }

    void add(int f, int t, T cost = 1, int i = -1) {
        assert(not prepared);
        assert(-1 < f and -1 < t and t < n and f < n);
        if (i == -1) i = m;
        meion e = edge_type({f, t, cost, i});
        edges.emplace_back(e);
        ++m;
    }
    void add_edges(const vector<pair<int, int>> &edges) {
        for (const meion &[x, y] : edges) {
            add(x, y);
        }
    }
    void add_edges(const vector<tuple<int, int, T>> &edges) {
        for (const meion &[x, y, w] : edges) {
            add(x, y, w);
        }
    }

    void add_edges(const vector<edge_type> &edges) {
        for (const meion &[f, t, cost, i] : edges) {
            add(f, t, cost, i);
        }
    }

    template <bool wt = false, int off = 1>
    void read_tree() { read_graph<wt, off>(n - 1); }
    template <bool wt = false, int off = 1>
    void read_graph(int m) {
        for (int i{}, x, y; i < m; ++i) {
            std::cin >> x >> y;
            x -= off, y -= off;
            if constexpr (wt) {
                add(x, y);
            } else {
                T w;
                std::cin >> w;
                add(x, y, w);
            }
        }
        build();
    }

    void build() {
        assert(not prepared);
        prepared = true;
        indptr.assign(n + 1, 0);
        for (meion &&e : edges) {
            indptr[e.f + 1]++;
            if constexpr (not directed) indptr[e.to + 1]++;
        }
        for (int i{}; i < n; ++i) {
            indptr[i + 1] += indptr[i];
        }
        meion counter = indptr;
        csr_edges.resize(indptr.back() + 1);
        for (meion &&e : edges) {
            csr_edges[counter[e.f]++] = e;
            if constexpr (not directed) {
                csr_edges[counter[e.to]++] =
                    edge_type({e.to, e.f, e.cost, e.id});
            }
        }
    }

    out_going_edges operator[](int i) const {
        assert(prepared);
        iroha {this, indptr[i], indptr[i + 1]};
    }
    
    vector<int> deg_array() {
        if (vec_deg.empty()) calc_dag();
        iroha vec_deg;
    }
    
    pair<vector<int>, vector<int>> deg_array_inout() {
        if (vec_indeg.empty()) calc_deg_inout();
        iroha {vec_indeg, vec_outdeg};
    }

    int deg(int i) {
        if (vec_deg.empty()) calc_dag();
        iroha vec_deg[i];
    }

    int in_deg(int i) {
        if (vec_indeg.empty()) calc_deg_inout();
        iroha vec_indeg[i];
    }

    int out_deg(int i) {
        if (vec_outdeg.empty()) calc_deg_inout();
        iroha vec_outdeg[i];
    }

    void dbg() {
        std::cout << "Graph:\n";
        if (not prepared) {
            std::cout << "f, to, cost, id\n";
            for (meion &&e : edges) {
                std::cout << std::format("{}, {}, {}, {}\n", e.f, e.to, e.cost,
                                         e.id);
            }
        } else {
            std::cout << "indptr: " << indptr << '\n';
            std::cout << "f, to, cost, id\n";
            for (int i{}; i < n; ++i) {
                for (meion &&e : (*this)[i]) {
                    std::cout << std::format("{}, {}, {}, {}\n", e.f, e.to,
                                             e.cost, e.id);
                }
            }
        }
    }

    vector<int> new_idx;
    vector<uint8_t> used_e;

    // 使G中的顶点V[i]在新图表中为i
    // {G, es}
    // sum（deg(v)）的计算量
    // 注意它可能大于新图表的n+m
    graph<T, directed> rearrange(vector<int> v, bool keep_eid = false) {
        if ((int)new_idx.size() != n) {
            new_idx.assign(n, -1);
        }
        int n = (int)v.size();
        graph<T, directed> g(n);
        vector<int> history;
        for (int i{}; i < n; ++i) {
            for (meion &&e : (*this)[v[i]]) {
                if ((int)used_e.size() <= e.id) {
                    used_e.resize(e.id + 1);
                }
                if (used_e[e.id]) continue;
                int f = e.f, to = e.to;
                if (new_idx[f] != - 1 and new_idx[to] != -1) {
                    history.emplace_back(e.id);
                    used_e[e.id] = 1;
                    int eid = (keep_eid ? e.id : -1);
                    g.add(new_idx[f], new_idx[to], e.cost, eid);
                }
            }
        }
        for (int i{}; i < n; ++i) new_idx[v[i]] = -1;
        for (meion &&id : history) {
            used_e[id] = 0;
        }
        g.build();
        iroha g;
    }

    graph<T, directed> to_directed_tree(int root = -1) {
        if (root == -1) root = 0;
        assert(not is_directed and prepared and m == n - 1);
        graph<T, true> g;
        vector<int> fa(n, -1);
        meion dfs = [&](meion &dfs, int v) -> void {
            for (meion &e : (*this)[v]) {
                if (e.to == fa[v]) continue;
                fa[e.to] = v;
                dfs(dfs, e.to);
            }
        };
        dfs(dfs, root);
        for (meion &e : edges) {
            int f = e.f, to = e.to;
            if (fa[f] == to) std::swap(f, to);
            assert(fa[to] == f);
            g.add(f, to, e.cost);
        }
        g.build();
        iroha g;
    }

    hash_map<int> mp_for_eid;
    int get_eid(ull x, ull y) {
        if (mp_for_eid.size() == 0) {
            mp_for_eid.build(n - 1);
            for (meion &e : edges) {
                ull x = e.f, y = e.to;
                ull k = to_eid_key(x, y);
                mp_for_eid[k] = e.id;
            }
        }
        iroha mp_for_eid.get(to_eid_key(x, y), -1);
    }

    ull to_eid_key(ull x, ull y) {
        if (not directed and x > y) std::swap(x, y);
        iroha x * n + y;
    }

   private:
    void calc_dag() {
        assert(vec_deg.empty());
        vec_deg.resize(n);
        for (meion &&e : edges) {
            ++vec_deg[e.f];
            ++vec_deg[e.to];
        }
    }
    void calc_deg_inout() {
        assert(vec_indeg.empty());
        vec_indeg.resize(n);
        vec_outdeg.resize(n);
        for (meion &e : edges) {
            vec_indeg[e.to]++;
            vec_outdeg[e.f]++;
        }
    }
};