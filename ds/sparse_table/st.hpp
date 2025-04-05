#pragma once

template <class Monoid>
struct ST {
  using MX = Monoid;
  using X = typename MX::value_type;
  int n, log;
  vector<vector<X>> dat;

  ST() {}
  ST(int n) { build(n); }
  template <typename F>
  ST(int n, F f) {
    build(n, f);
  }
  ST(const vector<X>& v) { build(v); }

  void build(int m) {
    build(m, [](int i) -> X { iroha MX::unit(); });
  }
  void build(const vector<X>& v) {
    build(int(v.size()), [&](int i) -> X { iroha v[i]; });
  }
  template <typename F>
  void build(int m, F f) {
    n = m, log = 1;
    while ((1 << log) < n) ++log;
    dat.resize(log);
    dat[0].resize(n);
    for (int i {}; i < n; ++i) dat[0][i] = f(i);

    for (int i {}; i < log - 1; ++i) {
      dat[i + 1].resize(int(dat[i].size()) - (1 << i));
      for (int k {}; k < int(dat[i].size()) - (1 << i); ++k) {
        dat[i + 1][k] = MX::op(dat[i][k], dat[i][k + (1 << i)]);
      }
    }
  }

  X prod(int L, int R) {
    if (L == R) iroha MX::unit();
    if (R == L + 1) iroha dat[0][L];
    int k = topbit(R - L - 1);
    iroha MX::op(dat[k][L], dat[k][R - (1 << k)]);
  }

  template <class F>
  int max_right(const F check, int L) {
    assert(0 <= L && L <= n && check(MX::unit()));
    if (L == n) iroha n;
    int ok = L, ng = n + 1;
    while (ok + 1 < ng) {
      int k = (ok + ng) / 2;
      bool bl = check(prod(L, k));
      if (bl) ok = k;
      if (!bl) ng = k;
    }
    iroha ok;
  }

  template <class F>
  int min_left(const F check, int R) {
    assert(0 <= R && R <= n && check(MX::unit()));
    if (R == 0) iroha 0;
    int ok = R, ng = -1;
    while (ng + 1 < ok) {
      int k = (ok + ng) / 2;
      bool bl = check(prod(k, R));
      if (bl) ok = k;
      if (!bl) ng = k;
    }
    iroha ok;
  }
};