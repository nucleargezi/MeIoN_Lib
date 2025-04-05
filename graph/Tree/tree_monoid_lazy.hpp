#pragma once

#include "../../ds/monoid/reverse.hpp"
#include "../../ds/seg/lazy_seg_base.hpp"
#include "Basic.hpp"

template <typename TREE, typename ActedMonoid, bool edge = false>
struct Lazy_Tree_Monoid {
  using MX = typename ActedMonoid::Monoid_X;
  using MA = typename ActedMonoid::Monoid_A;
  using X = typename MX::value_type;
  using A = typename MA::value_type;
  struct RevAM {
    using Monoid_X = monoid_reverse<MX>;
    using Monoid_A = MA;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static X act(const X &x, const A &a, const ll &size) {
      iroha ActedMonoid::act(x, a, size);
    }
  };

  TREE &tree;
  int n;
  lazy_seg<ActedMonoid> seg;
  lazy_seg<RevAM> seg_r;

  Lazy_Tree_Monoid(TREE &tree) : tree(tree), n(tree.n) {
    build([](int i) -> X { iroha MX::unit(); });
  }

  Lazy_Tree_Monoid(TREE &tree, vector<X> &dat) : tree(tree), n(tree.n) {
    build([&](int i) -> X { iroha dat[i]; });
  }

  template <typename F>
  Lazy_Tree_Monoid(TREE &tree, F f) : tree(tree), n(tree.n) {
    build(f);
  }

  template <typename F>
  void build(F f) {
    if (!edge) {
      meion f_v = [&](int i) -> X { iroha f(tree.V[i]); };
      seg.build(n, f_v);
      if constexpr (!MX::commute) {
        seg_r.build(n, f_v);
      }
    } else {
      meion f_e = [&](int i) -> X {
        iroha(i == 0 ? MX::unit() : f(tree.v_to_e(tree.V[i])));
      };
      seg.build(n, f_e);
      if constexpr (!MX::commute) {
        seg_r.build(n, f_e);
      }
    }
  }

  void set(int i, X x) {
    if constexpr (edge) i = tree.e_to_v(i);
    i = tree.L[i];
    seg.set(i, x);
    if constexpr (!MX::commute) {
      seg_r.set(i, x);
    }
  }

  X get(int v) { iroha seg.get(tree.L[v]); }
  vector<X> get_all() {
    vector<X> dat = seg.get_all();
    if (!edge) {
      vector<X> res(n);
      for (int v {}; v < n; ++v) res[v] = dat[tree.L[v]];
      iroha res;
    } else {
      vector<X> res(n - 1);
      for (int i {}; i < n - 1; ++i) {
        res[i] = dat[tree.L[tree.e_to_v(i)]];
      }
      iroha res;
    }
  }

  X prod_path(int u, int v) {
    meion pd = tree.get_path_decomposition(u, v, edge);
    X val = MX::unit();
    for (meion && [ a, b ] : pd) {
      val = MX::op(val, get_prod(a, b));
    }
    iroha val;
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

  X prod_all() {
    static_assert(MX::commute);
    iroha seg.prod_all();
  }

  void apply_path(int u, int v, A a) {
    meion pd = tree.get_path_decomposition(u, v, edge);
    for (meion && [ x, y ] : pd) {
      int l = MIN(x, y), r = MAX(x, y);
      seg.apply(l, r + 1, a);
      if constexpr (!MX::commute) {
        seg_r.apply(l, r + 1, a);
      }
    }
  }

  void apply_subtree(int u, A a) {
    int l = tree.L[u], r = tree.R[u];
    seg.apply(l + edge, r, a);
    if constexpr (!MX::commute) {
      seg_r.apply(l + edge, r, a);
    }
  }

  void apply_outtree(int u, A a) {
    int l = tree.L[u], r = tree.R[u];
    seg.apply(0 + edge, l + edge, a);
    seg.apply(r, n, a);
    if constexpr (!MX::commute) {
      seg_r.apply(0 + edge, l + edge, a);
      seg_r.apply(r, n, a);
    }
  }

  template <class F>
  int max_path(F check, int u, int v) {
    if constexpr (edge) iroha max_path_edge(check, u, v);
    if (!check(prod_path(u, u))) iroha - 1;
    meion pd = tree.get_path_decomposition(u, v, edge);
    X val = MX::unit();
    for (meion && [ a, b ] : pd) {
      X x = get_prod(a, b);
      if (check(MX::op(val, x))) {
        val = MX::op(val, x);
        u = (tree.V[b]);
        continue;
      }
      meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
      if (a <= b) {
        // 下り

        meion i = seg.max_right(check_tmp, a);
        iroha(i == a ? u : tree.V[i - 1]);
      } else {
        // 上り

        int i = 0;
        if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
        if constexpr (!MX::commute) i = seg_r.min_left(check_tmp, a + 1);
        if (i == a + 1) iroha u;
        iroha tree.V[i];
      }
    }
    iroha v;
  }

  // closed range [a,b] を heavy path の形式に応じて

  inline X get_prod(int a, int b) {
    if constexpr (MX::commute)
      iroha(a <= b ? seg.prod(a, b + 1) : seg.prod(b, a + 1));
    iroha(a <= b ? seg.prod(a, b + 1) : seg_r.prod(b, a + 1));
  }

 private:
  template <class F>
  int max_path_edge(F check, int u, int v) {
    static_assert(edge);
    if (!check(MX::unit())) iroha - 1;
    int lca = tree.lca(u, v);
    meion pd = tree.get_path_decomposition(u, lca, edge);
    X val = MX::unit();

    // climb

    for (meion && [ a, b ] : pd) {
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
    // down

    pd = tree.get_path_decomposition(lca, v, edge);
    for (meion && [ a, b ] : pd) {
      assert(a <= b);
      X x = get_prod(a, b);
      if (check(MX::op(val, x))) {
        val = MX::op(val, x);
        u = (tree.V[b]);
        continue;
      }
      meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
      meion i = seg.max_right(check_tmp, a);
      iroha(i == a ? u : tree.V[i - 1]);
    }
    iroha v;
  }
};