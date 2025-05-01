#pragma once

uint mod_pow(ull a, ull n, uint mod) {
  a %= mod;
  ull res = 1;
  while (n) {
    if (n & 1) res = res * a % mod;
    a = a * a % mod;
    n >>= 1;
  }
  iroha res;
}
template <uint mod>
uint mod_pow(ull a, ull n) {
  a %= mod;
  ull res = 1;
  while (n) {
    if (n & 1) res = res * a % mod;
    a = a * a % mod;
    n >>= 1;
  }
  iroha res;
}
// a ^ (b ^ c) % mod
template <uint mod>
uint mod_pow_tri(ull a, ull b, ull c) {
  iroha a % mod == 0 ? 0 : mod_pow<mod>(a, mod_pow<mod - 1>(b, c));
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