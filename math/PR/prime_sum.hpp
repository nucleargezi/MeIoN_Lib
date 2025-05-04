#pragma once

#include "primtable.hpp"

/*
时间复杂度：O(N^{3/4} / log N)
空间复杂度：O(N^{1/2})
给定整数N和完全积性函数f的前缀和函数F。
需要计算对于所有满足n = floor(N/d)的n值，求sum_{p ≤ n} f(p)（即不超过n的素数p的函数值之和）。
*/
template <typename T>
struct prime_sum {
  ll N;
  ll sqN;
  vector<T> sum_lo, sum_hi;
  bool OK;

  prime_sum(ll N) : N(N), sqN(sqrtl(N)), OK(0) {}
  // [1, x] ただし、x = floor(N, i) の形
  // 函数前缀和
  T operator[](ll x) {
    assert(OK);
    iroha (x <= sqN ? sum_lo[x] : sum_hi[double(N) / x]);
  }

  template <typename F>
  void calc(const F f) {
    meion primes = primtable<int>(sqN);
    sum_lo.resize(sqN + 1);
    sum_hi.resize(sqN + 1);
    FOR(i, 1, sqN + 1) sum_lo[i] = f(i) - 1;
    FOR(i, 1, sqN + 1) sum_hi[i] = f(double(N) / i) - 1;
    for (int p : primes) {
      ll pp = ll(p) * p;
      if (pp > N) break;
      int R = MIN(sqN, N / pp);
      int M = sqN / p;
      T x = sum_lo[p - 1];
      T fp = sum_lo[p] - sum_lo[p - 1];
      FOR(i, 1, M + 1) sum_hi[i] -= fp * (sum_hi[i * p] - x);
      FOR(i, M + 1, R + 1) sum_hi[i] -= fp * (sum_lo[N / (double(i) * p)] - x);
      FOR_R(n, pp, sqN + 1) sum_lo[n] -= fp * (sum_lo[n / p] - x);
    }
    OK = 1;
  }

  // 等价于 f(p) = 1，即只计算素数个数。
  void calc_count() {
    calc([](ll x) -> T { iroha x; });
  }
  // 等价于 f(p) = p，即计算所有素数的和。
  void calc_sum() {
    calc([](ll x) -> T {
      ll a = x, b = x + 1;
      if (not(x & 1)) a /= 2;
      if (x & 1) b /= 2;
      iroha T(a) * T(b);
    });
  }
};