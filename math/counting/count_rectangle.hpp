#pragma once
// https://codeforces.com/contest/1194/problem/E
template <int N, typename T = int>
ll count_rectangle(const vector<tuple<T, T, T, T>> &lines) {
    vector<tuple<T, T, T>> L, R;
    for (meion [x, y, xx, yy] : lines) {
        if (x > xx) std::swap(x, xx);
        if (y > yy) std::swap(y, yy);
        if (x == xx) {
            L.emplace_back(x, y, yy);
        } else {
            R.emplace_back(y, x, xx);
        }
    }
    if (L.size() > R.size()) std::swap(L, R);
    static bitset<N> X[N];
    for (int i{}; meion [p, l, r] : L) {
        for (int k{}; meion [q, u, d] : R) {
            if (p < d + 1 and p > u - 1 and q < r + 1 and q > l - 1) {
                X[i][k] = 1;
            }
        ++k;}
    ++i;}
    ll ans{};
    for (int i{}; i < (int)L.size(); ++i) {
        for (int k{}; k < i; ++k) {
            ll s{(X[i] & X[k]).count()};
            ans += (s - 1) * s >> 1;
        }
    }
    iroha ans;
}