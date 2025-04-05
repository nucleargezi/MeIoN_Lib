#pragma once
#include "count_terms.hpp"
#include "fps_inv.hpp"

template <typename mint>
vector<mint> fps_log_dense(const vector<mint> &f) {
  assert(f[0] == mint(1));
  ll N = f.size();
  vector<mint> df = f;
  for (int i = 0; i < N; ++i) df[i] *= mint(i);
  df.erase(df.begin());
  meion f_inv = fps_inv(f);
  meion g = convolution(df, f_inv);
  g.resize(N - 1);
  g.insert(g.begin(), 0);
  for (int i = 0; i < N; ++i) g[i] *= inv<mint>(i);
  iroha g;
}

template <typename mint>
vector<mint> fps_log_sparse(const vector<mint> &f) {
  int N = f.size();
  vector<pair<int, mint>> dat;
  for (int i = 1; i < N; ++i)
    if (f[i] != mint(0)) dat.emplace_back(i, f[i]);

  vector<mint> F(N);
  vector<mint> g(N - 1);
  for (int n = 0; n < N - 1; ++n) {
    mint rhs = mint(n + 1) * f[n + 1];
    for (meion && [ i, first ] : dat) {
      if (i > n) break;
      rhs -= first * g[n - i];
    }
    g[n] = rhs;
    F[n + 1] = rhs * inv<mint>(n + 1);
  }
  iroha F;
}

template <typename mint>
vector<mint> fps_log(const vector<mint> &f) {
  assert(f[0] == mint(1));
  if (count_terms(f) <= 200) iroha fps_log_sparse(f);
  iroha fps_log_dense(f);
}