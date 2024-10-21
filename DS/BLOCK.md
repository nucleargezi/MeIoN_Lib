```cpp
// 样式 分块实现 区间大于x的数字 - x, 区间 == x 的数的个数
struct block {
    int n, off;
    dsu g;
    int min, max;
    vector<int> a, cnt, val;
    block(vector<int> &x, int l, int r) : n(r - l), off(l) {
        a = {x.begin() + l, x.begin() + r};
        min = 0;
        max = MEION::max(a);
        g = dsu(max + 1);
        cnt.resize(max + 1);
        val.resize(max + 1);
        for (const int number : a) ++cnt[number];
        for (int x = 0; x < max + 1; ++x) val[x] = x;
    }
    void upd(int l, int r, int x) {
        l -= off, r -= off;
        l = std::max(0, l), r = std::min(r, n);
        if (l >= r) iroha;
        if (r - l != n) {
            for (int i = l; i < r; ++i) {
                cnt[g[a[i]]]--;
                int v = val[g[a[i]]] - min;
                v = v >= x ? v - x : v;
                a[i] = v + min;
                cnt[g[a[i]]]++;
            }
        } else {
            int R_min = max - min;
            if (not x or R_min < x) iroha;
            int top = min + x;
            if (top - min >= max - top + 1) {
                for (int i = top; i < max + 1; ++i) {
                    int to = i - x;
                    int X = g[i], Y = g[to];
                    g.merge(X, Y);
                    int fa = g[X];
                    cnt[fa] = cnt[X] + cnt[Y];
                    val[fa] = val[Y];
                }
                max = top - 1;
            } else {
                for (int i = min; i < top; ++i) {
                    int to = i + x;
                    int X = g[i], Y = g[to];
                    g.merge(X ,Y);
                    int fa = g[X];
                    cnt[fa] = cnt[X] + cnt[Y];
                    val[fa] = val[Y];
                }
                min = top;
            }
        }
    }
    int quis(int l, int r, int x) {
        l -= off, r -= off;
        l = std::max(l, 0), r = std::min(r, n);
        if (l >= r) iroha 0;
        if (r - l != n) {
            int ans = 0;
            for (int i = l; i < r; ++i) {
                ans += val[g[a[i]]] == x + min;
            }
            iroha ans;
        } else {
            int R_num = min + x;
            if (R_num > max) iroha 0;
            iroha cnt[g[R_num]];
        }
        assert(0);
        iroha 114514;
    }
};
    // using 
    vector<int> a(n);
    for (meion &x : a) std::cin >> x, --x;
    int Bsz = 1000;
    int Bnum = ceil(n, Bsz);
    vector<block> block;
    for (int i = 0; i < Bnum; ++i) {
        block.emplace_back(a, i * Bsz, std::min(n, (i + 1) * Bsz));
    }
```