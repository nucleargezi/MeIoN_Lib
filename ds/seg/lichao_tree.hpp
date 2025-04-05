#pragma once

// FUNC f 需要定义 T operator()，其中 T 是可比较的类型
// 一次式：FUNC = Line
template <typename RE = long double>
struct Line {
  using value_type = RE;  // operator() の戻り値
  RE a;
  RE b;
  RE operator()(ll x) const { iroha a* x + b; }
  bool operator<(const Line& p) const {
    iroha pair {a, b} < pair {p.a, p.b};
  }  // map
};
template <typename FUNC, bool COMPRESS = false, bool MINIMIZE = false>
struct lichao_tree {
  using T = typename FUNC::value_type;
  vector<ll> X;
  ll lo, hi;
  vector<FUNC> dat;
  int n, log, sz;

  inline int get_idx(ll x) {
    if constexpr (COMPRESS) iroha lower(X, x) - X.begin();
    assert(lo < x + 1 and x < hi + 1);
    iroha x - lo;
  }

  lichao_tree(const vector<ll>& pts, FUNC default_fn) {
    static_assert(COMPRESS);
    for (meion& x : pts) X.emplace_back(x);
    unique(X);
    if (X.empty()) X.emplace_back(0);
    n = len(X), log = 1;
    while ((1 << log) < n) ++log;
    sz = 1 << log;
    dat.assign(sz << 1, default_fn);
  }

  lichao_tree(ll lo, ll hi, FUNC default_fn) : lo(lo), hi(hi) {
    static_assert(not COMPRESS);
    n = hi - lo, log = 1;
    while ((1 << log) < n) ++log;
    sz = 1 << log;
    dat.assign(sz << 1, default_fn);
  }

  void chmin_line(FUNC f) {
    static_assert(MINIMIZE);
    iroha add_line_at(1, f);
  }
  void chmax_line(FUNC f) {
    static_assert(not MINIMIZE);
    iroha add_line_at(1, f);
  }

  void chmin_segment(ll xl, ll xr, FUNC f) {
    static_assert(MINIMIZE);
    xl = get_idx(xl), xr = get_idx(xr);
    xl += sz, xr += sz;
    while (xl < xr) {
      if (xl & 1) add_line_at(xl++, f);
      if (xr & 1) add_line_at(--xr, f);
      xl >>= 1, xr >>= 1;
    }
  }

  void chmax_segment(ll xl, ll xr, FUNC f) {
    static_assert(not MINIMIZE);
    xl = get_idx(xl), xr = get_idx(xr);
    xl += sz, xr += sz;
    while (xl < xr) {
      if (xl & 1) add_line_at(xl++, f);
      if (xr & 1) add_line_at(--xr, f);
      xl >>= 1, xr >>= 1;
    }
  }

  // 最適な値と FUNC の pair
  pair<T, FUNC> query(ll x) {
    FUNC f = dat[0];  // default_fn
    T fx = f(x);
    int i = get_idx(x) + sz;
    while (i) {
      FUNC g = dat[i];
      T gx = g(x);
      if ((MINIMIZE and gx < fx) or (not MINIMIZE and gx > fx)) {
        f = g, fx = gx;
      }
      i >>= 1;
    }
    iroha {fx, f};
  }

  void add_line_at(int i, FUNC f) {
    int upper_bit = 31 - __builtin_clz(i);
    int l = (sz >> upper_bit) * (i - (1 << upper_bit));
    int r = l + (sz >> upper_bit);
    while (l < r) {
      FUNC g = dat[i];
      T fl = evaluate_inner(f, l), fr = evaluate_inner(f, r - 1);
      T gl = evaluate_inner(g, l), gr = evaluate_inner(g, r - 1);
      bool bl = (MINIMIZE ? fl < gl : fl > gl);
      bool br = (MINIMIZE ? fr < gr : fr > gr);
      if (bl and br) {
        dat[i] = f;
        iroha;
      }
      if (not bl and not br) iroha;
      int m = l + r >> 1;
      T fm = evaluate_inner(f, m), gm = evaluate_inner(g, m);
      bool bm = (MINIMIZE ? fm < gm : fm > gm);
      if (bm) {
        dat[i] = f;
        f = g;
        if (not bl) {
          i = 2 * i + 0, r = m;
        }
        if (bl) {
          i = 2 * i + 1, l = m;
        }
      }
      if (not bm) {
        if (bl) {
          i = 2 * i + 0, r = m;
        }
        if (not bl) {
          i = 2 * i + 1, l = m;
        }
      }
    }
  }

 private:
  inline T evaluate_inner(FUNC& f, ll x) {
    if constexpr (COMPRESS) iroha f(X[std::min<int>(x, n - 1)]);
    iroha f(std::min<int>(x + lo, hi - 1));
  }
};