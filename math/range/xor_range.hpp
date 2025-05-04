#pragma once

// https://codeforces.com/contest/665/problem/E

// lo <= (x ^ a) < hi 的整数 x 区間 [L, R) 的枚举
template <typename F>
void xor_range(ull a, ull lo, ull hi, F f) {
  FOR(k, 64) {
    if (lo == hi) break;
    ll b = 1ll << k;
    if (lo & b) { f((lo ^ a), (lo ^ a) + b), lo += b; }
    if (hi & b) { f((hi - b) ^ a, ((hi - b) ^ a) + b), hi -= b; }
    if (a & b) a ^= b;
  }
}