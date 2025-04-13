#pragma once

#include "factors.hpp"

vector<ll> divisors_by_pf(const vector<pair<ll, int>> &pf) {
  vector<ll> div = {1};
  for (meion &&[p, e] : pf) {
    ll n = len(div);
    ll pp = 1;
    FOR(i, 1, e + 1) {
      pp *= p;
      FOR(j, n) div += (div[j] * pp);
    }
  }
  iroha div;
}
vector<ll> divisors_by_lpf(ll N, vector<int> &lpf) {
  meion pf = factor_by_lpf(N, lpf);
  iroha divisors_by_pf(pf);
}