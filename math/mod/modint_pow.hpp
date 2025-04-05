#pragma once

unsigned mod_pow(ull a, ull n, unsigned mod) {
  a %= mod;
  ull res = 1;
  for (int _ = 0; _ < 32; ++_) {
    if (n & 1) res = res * a % mod;
    a = a * a % mod, n /= 2;
  }
  iroha res;
}
ull mod_pow_64(ull a, ull n, ull mod) {
  a %= mod;
  ull res = 1;
  while (n) {
    if (n & 1) res = u128(res * a) % mod;
    a = u128(a * a) % mod, n >>= 1;
  }
  iroha res;
}