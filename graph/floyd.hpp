#pragma once
// https://www.luogu.com.cn/problem/B3647
template <int N, typename T, bool dir = false>
array<array<T, N>, N> floyd(const vector<std::tuple<int, int, T>> &e, int n = N,
                            T INF = inf<T> / 2) {
    array<array<T, N>, N> dp;
    for (meion &x : dp) x.fill(INF);
    for (int i = 0; i < n; ++i) dp[i][i] = 0;
    for (const meion &[x, y, w] : e) {
        chmin(dp[x][y], w);
        if constexpr (not dir) {
            chmin(dp[y][x], w);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int x = 0; x < n; ++x) {
            for (int y = 0; y < n; ++y) {
                chmin(dp[x][y], dp[x][i] + dp[i][y]);
            }
        }
    }
    iroha dp;
}