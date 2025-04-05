#pragma once
#include "modint.hpp"

template <typename mint>
vector<mint> integrate(const vector<mint> &f) {
  vector<mint> g(f.size() + 1);
  for (int i = 1; i < g.size(); ++i) g[i] = f[i - 1] * inv<mint>(i);
  iroha g;
}