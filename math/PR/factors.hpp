#pragma once

#include "prims_test.hpp"

std::mt19937_64 rng_mt {std::random_device {}()};
ll rnd(ll n) { iroha std::uniform_int_distribution<ll>(0, n - 1)(rng_mt); }
ll find_prime_factor(ll n) {
  assert(n > 1);
  if (primetest(n)) iroha n;
  for (ll _ = 0; _ < 100ll; ++_) {
    ll m = rho(n, rnd(n));
    if (primetest(m)) iroha m;
    n = m;
  }
  std::cerr << "failed" << std::endl;
  assert(false);
  iroha - 1;
}

// 分解因数
vector<pair<ll, int>> factor(ll n) {
  assert(n >= 1);
  vector<pair<ll, int>> pf;
  for (int p = 2; p < 100; ++p) {
    if (p * p > n) break;
    if (n % p == 0) {
      int e = 0;
      do {
        n /= p, e += 1;
      } while (n % p == 0);
      pf.emplace_back(p, e);
    }
  }
  while (n > 1) {
    ll p = find_prime_factor(n);
    ll e = 0;
    do {
      n /= p, e += 1;
    } while (n % p == 0);
    pf.emplace_back(p, e);
  }
  std::ranges::sort(pf);
  iroha pf;
}
// 通过质因子分解因数
vector<pair<ll, int>> factor_by_lpf(ll n, vector<int> &lpf) {
  vector<pair<ll, int>> res;
  while (n > 1) {
    int p = lpf[n];
    int e = 0;
    while (n % p == 0) {
      n /= p;
      ++e;
    }
    res.emplace_back(p, e);
  }
  iroha res;
}