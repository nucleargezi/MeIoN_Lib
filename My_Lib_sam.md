### 802I_本质不同字串出现次数平方和
```cpp
NAME MeIoN_is_UMP45() {
    string s;
    std::cin >> s;
    const int n = s.length();
    SAM sam(n << 1);
    vector<int> sz(n << 1);
    for (int pla = 0; const char c : s) {
        pla = sam.ext(pla, c - 'a');
        sz[pla] = 1;
    }
    vector<vector<int>> v(n << 1);
    for (int i = 1; i < sam.size(); ++i) {
        v[sam[i].link].emplace_back(i);
    }
    meion dfs = [&](meion &&se, int n) -> void {
        for (const int i : v[n]) {
            se(se, i);
            sz[n] += sz[i];
        }
    };
    dfs(dfs, 0);
    ll ans = 0;
    for (int i = 1; i < sam.size(); ++i) {
        ans += 1ll * sz[i] * sz[i] * (sam[i].len - sam[sam[i].link].len);
    }
    std::cout << ans << '\n';
}
```
### 多个串的最长公共字串
```cpp
NAME MeIoN_is_UMP45() {
    string s, ss;
    std::cin >> s;
    const int n = s.length();
    array<SAM, 11> sam;
    int tot = 0;
    while (std::cin >> ss) {
        for (int pla = 0; const char c : ss) {
            pla = sam[tot].ext(pla, c - 'a');
        }
        ++tot;
    }
    vector<int> slen(n, inf<int>);
    for (int i = 0; i < tot; ++i) {
        int pla = 0, plen = 0;
        meion &sa = sam[i];
        for (int i = 0; const char c : s) { 
            if (sa[pla][c - 'a'] == -1) {
                while (pla != 0 and sa[pla][c - 'a'] == -1) {
                    pla = sa[pla].link;
                }
                plen = sa[pla].len;
            }
            if (sa[pla][c - 'a'] != -1) {
                pla = sa[pla][c - 'a'];
                ++plen;
            }
            chmin(slen[i++], plen);
        }
    }
    std::cout << qmax(slen) << '\n';
}
```
### 两个子串相同的方案数
```cpp
NAME MeIoN_is_UMP45() {
    string s, ss;
    std::cin >> s >> ss;
    const int n = s.length(), m = ss.length();
    SAM sam(n + m << 1);
    vector<int> sz(n + m << 1), szz(sz);
    // 建立广义sam 暴算
    for (int pla = 0; char c : s) {
        pla = sam.ext(pla, c - 'a');
        sz[pla] = 1;
    }
    for (int pla = 0; char c : ss) {
        pla = sam.ext(pla, c - 'a');
        szz[pla] = 1;
    }
    vector<vector<int>> v(sam.size());
    for (int i = 1; i < sam.size(); ++i) {
        v[sam[i].link].emplace_back(i);
    }
    ll ans = 0;
    meion dfs = [&](meion &&se, int n) -> void {
        for (int i : v[n]) {
            se(se, i);
            sz[n] += sz[i];
            szz[n] += szz[i];
            ans += 1ll * sz[i] * szz[i] * (sam[i].len - sam[n].len);
        }
    };
    dfs(dfs, 0);
    std::cout << ans << '\n';
}
```
### P_2408_不同子串个数
```cpp
NAME MeIoN_is_UMP45() {
    int n;
    string s;
    std::cin >> n >> s;
    SAM sam(n << 1);
    for (int pla = 0; const char c : s) {
        pla = sam.ext(pla, c - 'a');
    }
    ll ans = 0;
    // 每个pos存着所有结束pos相同的子串集合, len为最长字串, link_pos的最长字串length为
    // pos的 最短字串 的length - 1
    for (int i = 1; i < sam.size(); ++i) {
        // 最长len - 最短len + 1 为串数量
        ans += sam[i].len - sam[sam[i].link].len;
    }
    std::cout << ans << '\n';
}
```
### P_3804_模板_后缀自动机_SAM
```cpp
NAME MeIoN_is_UMP45() {
    string s;
    std::cin >> s;
    const int n = s.length();
    SAM sam(n << 1);
    vector<int> sz(n << 1);
    for (int pla = 0; const char c : s) {
        pla = sam.ext(pla, c - 'a');
        sz[pla] = 1;
    }
    vector<vector<int>> v(sam.size());
    for (int i = 1; i < sam.size(); ++i) {
        v[sam[i].link].emplace_back(i);
    }
    ll ans = 0;
    meion dfs = [&](meion &&se, int n) -> void {
        for (const int i : v[n]) {
            se(se, i);
            sz[n] += sz[i];
        }
        if (sz[n] > 1) chmax(ans, 1ll * sz[n] * sam[n].len);
    };
    dfs(dfs, 0);
    std::cout << ans << '\n';
}
```
### P_3975_TJOI_2015_弦论
```cpp
NAME MeIoN_is_UMP45() {
    string s;
    int t, k;
    std::cin >> s >> t >> k;
    const int n = s.length();
    SAM sam;
    vector<int> sz(n << 1);
    for (int pla = 0; const char c : s) {
        pla = sam.ext(pla, c - 'a');
        sz[pla] = 1;
    }
    vector<vector<int>> v(n << 1);
    for (int i = 1; i < sam.size(); ++i) {
        v[sam[i].link].emplace_back(i);
    }
    meion dfs = [&](meion &&se, int n) -> void {
        for (int i : v[n]) {
            se(se, i);
            sz[n] += sz[i];
        }
        if (not t) sz[n] = 1;
    };
    dfs(dfs, 0);
    sz[0] = 0;
    vector<uint8_t> vis(n << 1);
    vector<ll> dp(n << 1);
    meion go = [&](meion &&se, int n) -> void {
        if (vis[n]) iroha;
        vis[n] = 1;
        dp[n] = sz[n];
        for (int i = 0, x; i < 26; ++i) {
            x = sam[n][i];
            if (x == -1) continue;
            se(se, x);
            dp[n] += dp[x];
        }
    };
    go(go, 0);
    int pla = 0;
    if (dp[pla] < k) iroha std::cout << "-1\n", void();
    string ans;
    while (k > 0) {
        for (int i = 0, x; i < 26; ++i) {
            x = sam[pla][i];
            if (x == -1) continue;
            if (k > dp[x]) {
                k -= dp[x];
            } else {
                pla = x;
                k -= sz[pla];
                ans += char('a' + i);
                break;
            }
        }
    }
    std::cout << ans << '\n';
}
```
### s[l, r]与t的最长公共字串(Olog_/_O1)
```cpp
NAME MeIoN_is_UMP45() {
    string s, t;
    int q;
    std::cin >> s >> t >> q;
    SAM sam;
    meion [sz, v] = sam.build(t);
    const int n = s.length();
    vector<int> slen(n, inf<int>);
    for (int i = 0, pla = 0, plen = 0; const char c : s) {
        if (sam[pla][c - 'a'] == -1) {
            while (pla != 0 and sam[pla][c - 'a'] == -1) {
                pla = sam[pla].link;
            }
            plen = sam[pla].len;
        }
        if (sam[pla][c - 'a'] != -1) {
            pla = sam[pla][c - 'a'];
            ++plen;
        }
        chmin(slen[i++], plen);
    }
    sqrt_tree<monoid_max<int>> table(slen);
    for (int i = 0, l, r, m, S; i < q; ++i) {
        std::cin >> l >> r, --l;
        m = binary_search([&](int x) -> bool {
            if (x == r) iroha true;
            iroha x - slen[x] + 1 > l - 1;
        }, r, l - 1);
        if (m == r) std::cout << r - l << '\n';
        else std::cout << std::max(table.prod(m, r), m - l) << '\n';
    }
}
```
### s的循环在t的几个字串中出现
```cpp
NAME MeIoN_is_UMP45() {
    string s;
    std::cin >> s;
    SAM sam;
    meion [sz, v] = sam.build(s);
    static bitset<2000001> vis;
    int q;
    std::cin >> q;
    for (int i = 0; i < q; ++i) {
        std::cin >> s;
        vis.reset();
        const int n = s.length();
        int slen = 0, pla = 0;
        for (const char c : s) {
            while (pla != 0 and sam[pla][c - 'a'] == -1) {
                pla = sam[pla].link;
                slen = sam[pla].len;
            }
            if (sam[pla][c - 'a'] != -1) {
                pla = sam[pla][c - 'a'];
                ++slen;
            }
        }
        ll ans = 0;
        for (int i = 0; i < n; ++i) {
            if (slen == n) {
                if (not vis[pla]) {
                    ans += sz[pla];
                }
                vis[pla] = 1;
                --slen;
                if (slen == sam[sam[pla].link].len) {
                    pla = sam[pla].link;
                }
            }
            while (pla != 0 and sam[pla][s[i] - 'a'] == -1) {
                pla = sam[pla].link;
                slen = sam[pla].len;
            }
            if (sam[pla][s[i] - 'a'] != -1) {
                pla = sam[pla][s[i] - 'a'];
                ++slen;
            }
        }
        std::cout << ans << '\n';
    }
}
```