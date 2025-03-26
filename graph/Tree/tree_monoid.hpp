#pragma once
#include "../../ds/seg/seg_base.hpp"
#include "../../ds/monoid/reverse.hpp"
#include "Basic.hpp"

// P4427 [BJOI2018] 求和 prod_path

template <typename TREE, typename monoid, bool edge = false>
struct tree_monoid {
    using MX = monoid;
    using X = typename MX::value_type;
    TREE &tree;
    int n;
    Seg<MX> seg;
    Seg<monoid_reverse<MX>> seg_r;

    tree_monoid(TREE &tree) : tree(tree), n(tree.n) {
        build([](int i) -> X { iroha MX::unit(); });
    }
    tree_monoid(TREE &tree, vector<X> &dat) : tree(tree), n(tree.n) {
        build([&](int i) -> X { iroha dat[i]; });
    }
    template <typename F>
    tree_monoid(TREE &tree, F f) : tree(tree), n(tree.n) {
        build(f);
    }
    template <typename F>
    void build(F f) {
        if (not edge) {
            meion f_v = [&](int i) -> X { iroha f(tree.V[i]); };
            seg.build(n, f_v);
            if constexpr (not MX::commute) {
                seg_r.build(n, f_v);
            }
        } else {
            meion f_e = [&](int i) -> X {
                iroha (i == 0 ? MX::unit() : f(tree.v_to_e(tree.V[i])));
            };
            seg.build(n, f_e);
            if constexpr (not MX::commute) {
                seg_r.build(n, f_e);
            }
        }
    }

    void set(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.set(i, x);
        if constexpr (not MX::commute) seg_r.set(i, x);
    }

    void multiply(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.multiply(i, x);
        if constexpr (not MX::commute) seg_r.multiply(i, x);
    }
    void apply(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.multiply(i, x);
        if constexpr (not MX::commute) seg_r.multiply(i, x);
    }

    X prod_path(int u, int v) {
        meion pd = tree.get_path_decomposition(u, v, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            val = MX::op(val, get_prod(a, b));
        }
        iroha val;
    }

    // 在 uv 路径上，找到满足 check 条件的 prod_path(u, x) 的最后一个 x。
    // 如果没有找到（即 path(u, u) 不符合要求），返回 -1。
    template <class F>
    int max_path(F check, int u, int v) {
        if constexpr (edge) iroha max_path_edge(check, u, v);
        if (not check(prod_path(u, u))) iroha -1;
        meion pd = tree.get_path_decomposition(u, v, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.V[b]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            if (a <= b) {
                meion i = seg.max_right(check_tmp, a);
                iroha (i == a ? u : tree.V[i - 1]);
            } else {
                int i = 0;
                if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
                if constexpr (not MX::commute)
                    i = seg_r.min_left(check_tmp, a + 1);
                if (i == a + 1) iroha u;
                iroha tree.V[i];
            }
        }
        iroha v;
    }

    X prod_subtree(int u, int root = -1) {
        if (root == u) iroha prod_all();
        if (root == -1 || tree.in_subtree(u, root)) {
            int l = tree.L[u], r = tree.R[u];
            iroha seg.prod(l + edge, r);
        }
        assert(!edge);  // さぼり

        u = tree.jump(u, root, 1);
        int L = tree.L[u], R = tree.R[u];
        iroha MX::op(seg.prod(0, L), seg.prod(R, n));
    }

    X prod_all() { iroha prod_subtree(tree.V[0]); }

    inline X get_prod(int a, int b) {
        if constexpr (MX::commute) {
            iroha (a <= b) ? seg.prod(a, b + 1) : seg.prod(b, a + 1);
        }
        iroha (a <= b) ? seg.prod(a, b + 1) : seg_r.prod(b, a + 1);
    }

   private:
    template <class F>
    int max_path_edge(F check, int u, int v) {
        static_assert(edge);
        if (!check(MX::unit())) iroha -1;
        int lca = tree.lca(u, v);
        meion pd = tree.get_path_decomposition(u, lca, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            assert(a >= b);
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.parent[tree.V[b]]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            int i = 0;
            if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
            if constexpr (!MX::commute) i = seg_r.min_left(check_tmp, a + 1);
            if (i == a + 1) iroha u;
            iroha tree.parent[tree.V[i]];
        }
        pd = tree.get_path_decomposition(lca, v, edge);
        for (meion &&[a, b] : pd) {
            assert(a <= b);
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.V[b]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            meion i = seg.max_right(check_tmp, a);
            iroha (i == a ? u : tree.V[i - 1]);
        }
        iroha v;
    }
};