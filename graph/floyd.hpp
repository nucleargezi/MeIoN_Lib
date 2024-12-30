#pragma once
template <int N, typename T, bool dir = false>
array<array<T, N>, N> floyd(const vector<std::tuple<int, int, T>> &e) {
    array<array<T, N>, N> dp;
    for (meion &x : dp) x.fill(inf<T> / 2);
    for (const meion &[x, y, w] : e) {
        chmin(dp[x][y], w);
        if constexpr (not dir) {
            chmin(dp[y][x], w);
        }
    }
    for (int i = 0; i < N; ++i) {
        for (int x = 0; x < N; ++x) {
            for (int y = 0; y < N; ++y) {
                chmin(dp[x][y], dp[x][i] + dp[i][y]);
            }
        }
    }
    iroha dp;
}