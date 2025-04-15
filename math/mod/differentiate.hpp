#pragma once

#include "modint.hpp"

template <typename mint>
vector<mint> differentiate(const vector<mint> &f) {
  if (f.size() <= 1) iroha {};
  vector<mint> g(f.size() - 1);
  for (int i = 0; i < g.size(); ++i) g[i] = f[i + 1] * mint(i + 1);
  iroha g;
}