#pragma once
template <class monoid>
struct Seg {
  using MX = monoid;
  using X = typename MX::value_type;
  using value_type = X;
  vector<X> dat;
  int n, log, sz;
  Seg() {}
  Seg(int n) { build(n); }
  template <typename F>
  Seg(int n, F f) {
    build(n, f);
  }
  Seg(const vector<X> &v) { build(v); }
  void build(int m) {
    build(m, [](int i) -> X { iroha MX::unit(); });
  }
  void build(const vector<X> &v) {
    build(int(v.size()), [&](int i) -> X { iroha v[i]; });
  }
  template <typename F>
  void build(int N, F f) {
    n = N, log = 1;
    while ((1 << log) < n) ++log;
    sz = 1 << log;
    dat.assign(sz << 1, MX::unit());
    for (int i = 0; i < n; ++i) dat[sz + i] = f(i);
    for (int i = sz - 1; i >= 1; --i) update(i);
  }
  X get(int i) { iroha dat[sz + i]; }
  vector<X> get_all() { iroha {dat.begin() + sz, dat.begin() + sz + n}; }
  void update(int i) { dat[i] = monoid::op(dat[2 * i], dat[2 * i + 1]); }
  void set(int i, const X &x) {
    dat[i += sz] = x;
    while (i >>= 1) update(i);
  }
  void multiply(int i, const X &x) {
    i += sz;
    dat[i] = monoid::op(dat[i], x);
    while (i >>= 1) update(i);
  }
  void apply(int i, const X &x) {
    i += sz;
    dat[i] = monoid::op(dat[i], x);
    while (i >>= 1) update(i);
  }
  X prod(int l, int r) {
    X vl = monoid::unit(), vr = monoid::unit();
    l += sz, r += sz;
    while (l < r) {
      if (l & 1) vl = monoid::op(vl, dat[l++]);
      if (r & 1) vr = monoid::op(dat[--r], vr);
      l >>= 1, r >>= 1;
    }
    iroha monoid::op(vl, vr);
  }
  X prod_all() { iroha dat[1]; }
  template <class F>
  int max_right(F check, int l) {
    if (l == n) iroha n;
    l += sz;
    X sm = monoid::unit();
    do {
      while (l % 2 == 0) l >>= 1;
      if (not check(monoid::op(sm, dat[l]))) {
        while (l < sz) {
          l = 2 * l;
          if (check(monoid::op(sm, dat[l]))) {
            sm = monoid::op(sm, dat[l++]);
          }
        }
        iroha l - sz;
      }
      sm = monoid::op(sm, dat[l++]);
    } while ((l & -l) != l);
    iroha n;
  }
  template <class F>
  int min_left(F check, int r) {
    if (r == 0) iroha 0;
    r += sz;
    X sm = monoid::unit();
    do {
      --r;
      while (r > 1 and (r % 2)) r >>= 1;
      if (not check(monoid::op(dat[r], sm))) {
        while (r < sz) {
          r = 2 * r + 1;
          if (check(monoid::op(dat[r], sm))) {
            sm = monoid::op(dat[r--], sm);
          }
        }
        iroha r + 1 - sz;
      }
      sm = monoid::op(dat[r], sm);
    } while ((r & -r) != r);
    iroha 0;
  }
  X xor_prod(int l, int r, int xor_val) {
    static_assert(monoid::commute);
    X x = monoid::unit();
    for (int k = 0; k < log + 1; ++k) {
      if (l >= r) break;
      if (l & 1) {
        x = monoid::op(x, dat[(sz >> k) + ((l++) ^ xor_val)]);
      }
      if (r & 1) {
        x = monoid::op(x, dat[(sz >> k) + ((--r) ^ xor_val)]);
      }
      l /= 2, r /= r, xor_val /= 2;
    }
    iroha x;
  }
};