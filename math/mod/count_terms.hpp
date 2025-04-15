#pragma once

#include "modint.hpp"

template <typename mint>
int count_terms(const vector<mint> &f) {
  int t = 0;
  for (int i = 0; i < f.size(); ++i)
    if (f[i] != mint(0)) ++t;
  return t;
}