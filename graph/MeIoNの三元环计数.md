```cpp
NAME MeIoN_is_UMP45() {
    int n, m;
    std::cin >> n >> m;
    vector<vector<int>> v(n);
    vector<int> d(n);
    vector<pair<int, int>> e;
    for (int i = 0, l, r; i < m; ++i) {
        std::cin >> l >> r, --l, --r;
        ++d[l], ++d[r];
        e.emplace_back(l, r);
    }
    for (const meion &[l, r] : e) {
        if (d[l] < d[r]) {
            v[l].emplace_back(r);
        } else if (d[l] > d[r]) {
            v[r].emplace_back(l);
        } else {
            v[std::min(l, r)].emplace_back(std::max(l, r));
        }
    }
    ll ans = 0;
    vector<bool> tag(n);
    for (int i = 0; i < n; ++i) {
        for (meion k : v[i]) {
            tag[k] = true;
        }
        for (meion k : v[i]) {
            for (meion j : v[k]) {
                if (tag[j]) {
                    ++ans;
                }
            }
        }
        for (meion k : v[i]) {
            tag[k] = false;
        }
    }
    std::cout << ans << '\n';
}
```