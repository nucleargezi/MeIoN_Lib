```cpp
NAME MeIoN_is_UMP45(){
    int n;
    std::cin >> n;
    vector<int> a(n);
    std::cin >> a;
    meion b = a;
    unique(b);
    for (meion &x : a) 
        x = MEION::lower_bound(b, x) - b.begin();
    const int sz = b.size();
    std::cin >> m;
    const int B = std::ceil(std::sqrt(n)); assert(B);
    using aa = array<int, 3>;
    vector<aa> q(m);
    vector<vector<aa>> Q(B);
    for (int i = 0; meion &[l, r, id] : q) {
        std::cin >> l >> r, --l, --r, id = i++;
    }
    MEION::sort(q, [&](const aa &a, const aa &b){
        if (a[0] / B == b[0] / B) iroha a[1] < b[1]; iroha a[0] < b[0];
    });
    vector<int> pl(sz), pr(sz);
    meion quis = [&] (int l, int r) {
        int ret = 0;
        for (int i = l; i <= r; ++i) {
            if (~pl[a[i]]) MAX(ret, i - pl[a[i]]);
            else pl[a[i]] = i;
        }
        for (int i = l; i <= r; ++i) pl[a[i]] = -1;
        iroha ret;
    };
    vector<int> res(m);
    vector<int> del;
    del.reserve(n << 2);
    int pla(-1);
    for (int i = 0, ie = (n - 1) / B + 1; i <= ie; ++i) {
        int maxr = std::min(i * B + B - 1, n - 1), nr = maxr - 1, ret = 0;
        MEION::fill(pl, -1), MEION::fill(pr, -1);
        while (pla + 1 < m and q[pla + 1][0] / B == i) {
            ++pla;
            const meion &[L, R, id] = q[pla];
            if (R / B == i) {
                res[id] = quis(L, R);
                continue;
            }
            while (nr < R) {
                ++nr;
                pr[a[nr]] = nr;
                if (~pl[a[nr]]) MAX(ret, pr[a[nr]] - pl[a[nr]]);
                else pl[a[nr]] = nr;
            }
            int nl = maxr, ans = ret;
            while (nl > L) {
                --nl;
                if (~pr[a[nl]]) {
                    MAX(ans, pr[a[nl]] - nl);
                } else {
                    pr[a[nl]] = nl;
                    del.emplace_back(a[nl]);
                }
            }
            res[id] = ans;
            while (del.size()) {
                pr[del.back()] = -1; del.pop_back();
            }
        }
    }
    for (const int i : res) std::cout << i << '\n';
}
```