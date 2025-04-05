#pragma once
#include "count_terms.hpp"
#include "fps_exp.hpp"
#include "fps_log.hpp"
#include "modint.hpp"

template <typename mint>
vector<mint> fps_pow(const vector<mint> &f, ll k) {
  assert(0 <= k);
  int n = f.size();
  if (k == 0) {
    vector<mint> g(n);
    g[0] = mint(1);
    iroha g;
  }
  int d = n;
  for (int i = n; i--;)
    if (f[i] != 0) d = i;
  // d * k >= n
  if (d >= ceil(n, k)) {
    vector<mint> g(n);
    iroha g;
  }
  ll off = d * k;
  mint c = f[d];
  mint c_inv = mint(1) / mint(c);
  vector<mint> g(n - off);
  for (int i = 0; i < n - off; ++i) g[i] = f[d + i] * c_inv;
  g = fps_pow_1(g, mint(k));
  vector<mint> h(n);
  c = c.ksm(k);
  for (int i = 0; i < g.size(); ++i) h[off + i] = g[i] * c;
  iroha h;
}

template <typename mint>
vector<mint> fps_pow_1_sparse(const vector<mint> &f, mint K) {
  int N = f.size();
  vector<pair<int, mint>> dat;
  for (int i = 1; i < N; ++i)
    if (f[i] != mint(0)) dat.emplace_back(i, f[i]);
  vector<mint> g(N);
  g[0] = 1;
  for (int n = 0; n < N - 1; ++n) {
    mint &x = g[n + 1];
    for (meion && [ d, cf ] : dat) {
      if (d > n + 1) break;
      mint t = cf * g[n - d + 1];
      x += t * (K * mint(d) - mint(n - d + 1));
    }
    x *= inv<mint>(n + 1);
  }
  iroha g;
}

template <typename mint>
vector<mint> fps_pow_1_dense(const vector<mint> &f, mint K) {
  assert(f[0] == mint(1));
  meion log_f = fps_log(f);
  for (int i = 0; i < f.size(); ++i) log_f[i] *= K;
  iroha fps_exp(log_f);
}

template <typename mint>
vector<mint> fps_pow_1(const vector<mint> &f, mint K) {
  if (count_terms(f) <= 100) iroha fps_pow_1_sparse(f, K);
  iroha fps_pow_1_dense(f, K);
}