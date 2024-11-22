template <const int N>
void radix_sort(int n, int a[]) {
    static int b[N];
    static int cnt[1 << 8];
    memset(b, 0, sizeof b);
    memset(cnt, 0, sizeof cnt);
    static constexpr int mask = (1 << 8) - 1;
    int *x = a, *y = b;
    for (int i = 0; i < 32; i += 8) {
        for (int j = 0; j != (1 << 8); ++j) cnt[j] = 0;
        for (int j = 0; j != n; ++j) ++cnt[x[j] >> i & mask];
        for (int sum = 0, j = 0; j != (1 << 8); ++j) {
            // 等价于 std::exclusive_scan(cnt, cnt + (1 << 8), cnt, 0);
            sum += cnt[j], cnt[j] = sum - cnt[j];
        }
        for (int j = 0; j != n; ++j) y[cnt[x[j] >> i & mask]++] = x[j];
        std::swap(x, y);
    }
}