#pragma once
template <int N, typename T, bool dir = false>
array<array<T, N>, N> floyd(const vector<std::tuple<int, int, T>> &e, int n = N) {
    array<array<T, N>, N> dp;
    for (meion &x : dp) x.fill(inf<T> / 2);
    for (int i = 0; i < n; ++i) dp[i][i] = 0;
    for (const meion &[x, y, w] : e) {
        chmin(dp[x][y], w);
        if constexpr (not dir) {
            chmin(dp[y][x], w);
        }
    }
    for (int i = 0; i < n + 1; ++i) {
        for (int x = 0; x < n + 1; ++x) {
            for (int y = 0; y < n + 1; ++y) {
                chmin(dp[x][y], dp[x][i] + dp[i][y]);
            }
        }
    }
    iroha dp;
}