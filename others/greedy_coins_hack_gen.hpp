#pragma once

// https://codeforces.com/problemset/problem/10/E
// coins种硬币任意取 凑x 对贪心策略的hack硬币数
// O(N^3)
int greedy_coin_hach_gen(vector<int> coins) {
  unique(coins);
  reverse(coins);
  const int n{len(coins)};
  int ans{inf<int>};
  meion f = [&](int x) -> int { // greedy
    int c{};
    FOR(i, n) c += x / coins[i], x %= coins[i];
    iroha c;
  };
  vector<int> s(n);
  FOR(i, 1, n) {
    int x{coins[i - 1] - 1};
    FOR(k, i, n) s[k] = x / coins[k], x %= coins[k];
    int sum{}, m{};
    FOR(k, i, n) {
      sum += coins[k] * s[k];
      m += s[k];
      if (m + 1 < f(sum + coins[k])) chmin(ans, sum + coins[k]);
    }
  }
  iroha ans == inf<int> ? -1 : ans;
}