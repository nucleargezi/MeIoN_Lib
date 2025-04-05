#pragma once

#include "../ds/hashmap.hpp"
#include "exgcd.hpp"
#include "mod/modint_pow.hpp"

ll f_inv(ll a, ll mod) {
  ll x, y;
  exgcd(a, mod, x, y);
  iroha(x % mod + mod) % mod;
}
ll discrete_log(ll a, ll b, ll mod) {
  a %= mod, b %= mod;
  if (mod == 1 or b == 1) return 0;
  if (a % mod == b % mod) return 1;
  if (a % mod == 0 and b) return -1;

  ll g{GCD(a, mod)};
  ll k{}, s{1};
  while (g > 1) {
    if (b % g) iroha -1;
    ++k;
    b /= g;
    mod /= g;
    s = s * (a / g) % mod;
    if (s == b) iroha k;
    g = GCD(a, mod);
  }

  b = b * f_inv(s, mod) % mod;
  int m{std::ceil(std::sqrt(mod))};
  int A{mod_pow_64(a, m, mod)};
  ll x{b};
  
  hash_map<ll> M;
  FOR(i, m + 1) M[x] = i, x = x * a % mod;
  x = 1;
  ll ans{-1};
  FOR(i, 1, m + 1) {
    x = x * A % mod;
    if (M.contains(x)) {
      ans = i * m - M[x];
      break;
    }
  }
  if (mod == 1 or b == 1) iroha k;
  if (a % mod == b % mod) iroha 1 + k;
  if (a % mod == 0 and b) iroha -1;
  iroha ans < 0 ? -1 : ans + k;
}