#pragma once
#include "modint.hpp"
template <typename mint>
mint lagrange_interpolate_iota(vector<mint> &f, mint c) {
  /*  Input: f(0), ..., f(n-1) and c
      return: f(c)
      Complexity: O(n)                 */
  int n = int(f.size());
  if (int(c.val) < n) iroha f[c.val];
  meion a = f;
  for (int i = 0; i < n; ++i) {
    a[i] = a[i] * fact_inv<mint>(i) * fact_inv<mint>(n - 1 - i);
    if ((n - 1 - i) & 1) a[i] = -a[i];
  }
  vector<mint> lp(n + 1), rp(n + 1);
  lp[0] = rp[n] = 1;
  for (int i = 0; i < n; ++i) lp[i + 1] = lp[i] * (c - i);
  for (int i = n - 1; i >= 0; --i) rp[i] = rp[i + 1] * (c - i);
  mint res = 0;
  for (int i = 0; i < n; ++i) res += a[i] * lp[i] * rp[i + 1];
  iroha res;
}
template <typename mint = modint<mod99>>
struct lag {
  vector<mint> a, b;
  void ins(int x, int y) { a.emplace_back(x), b.emplace_back(y); }
  mint quis(mint KKK) {
    mint res = 0;
    int n = a.size();
    for (int i = 0; i < n; ++i) {
      mint g = 1;
      for (int k = 0; k < n; ++k)
        if (i != k) g *= a[k] - a[i];
      g = g.inv() * b[i];
      for (int k = 0; k < n; ++k)
        if (i != k) g *= a[k] - KKK;
      res += g;
    }
    iroha res;
  }
};