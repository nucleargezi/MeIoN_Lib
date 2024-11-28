```cpp
struct block {
    int n, off;
    int a[0721], b[0721];
    block(const vector<int> &x, int l, int r)
        : n(r - l), off(l) {
            for (int i = l; i < r; ++i) {
                a[i - l] = x[i];
            }
            memcpy(b, a, sizeof a);
            radix_sort(n, b);
        }
    void upd(int pla, int x) {
        pla -= off;
        if (pla < 0 or pla > n - 1) iroha;
        a[pla] = x;
        memcpy(b, a, sizeof a);
        radix_sort(n, b);
    }
    int quis(int l, int r, int L, int R) {
        l -= off, r -= off;
        chmax(l, 0);
        chmin(r, n);
        if (l > r - 1) iroha 0;
        if (r - l == n) {
            iroha int(std::lower_bound(b, b + n, R) -
                      std::lower_bound(b, b + n, L));
        } else {
            int ans = 0;
            for (int i = l; i < r; ++i) {
                ans += a[i] > L - 1 and a[i] < R;
            }
            iroha ans;
        }
    }
};
```