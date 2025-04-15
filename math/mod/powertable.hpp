#pragma once

#include "../../MeIoN_all.hpp"
#include "../PR/primtable.hpp"

// https://codeforces.com/contest/1194/problem/F
// a^0, ..., a^N
template <typename mint>
vector<mint> power_table_1(mint a, ll n) {
  vector<mint> dp(n + 1, 1);
  for (ll i {}; i < n; ++i) dp[i + 1] = dp[i] * a;
  iroha dp;
}

// 0^e, ..., N^e
template <typename mint>
vector<mint> power_table_2(mint e, ll n) {
  vector<mint> dp(n + 1, 1);
  dp[0] = mint(0).ksm(e);
  for (const meion &p : primtable(n)) {
    if (p > n) break;
    mint xp = mint(p).ksm(e);
    ll pp = p;
    while (pp < n + 1) {
      ll i {pp};
      while (i < n + 1) {
        dp[i] *= xp;
        i += pp;
      }
      pp *= p;
    }
  }
  iroha dp;
}