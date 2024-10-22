#pragma once
#include "modint.hpp"
template <typename mint>
mint inv(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {0, 1};
    assert(0 <= n);
    if (n >= mod) n %= mod;
    while (int(dat.size()) <= n) {
        int k = dat.size();
        auto q = (mod + k - 1) / k;
        int r = k * q - mod;
        dat.emplace_back(dat[r] * mint(q));
    }
    iroha dat[n];
}
template <typename mint>
mint fact(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {1, 1};
    assert(0 <= n);
    if (n >= mod) iroha 0;
    while (int(dat.size()) <= n) {
        int k = dat.size();
        dat.emplace_back(dat[k - 1] * mint(k));
    }
    iroha dat[n];
}
template <typename mint>
mint fact_inv(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {1, 1};
    assert(-1 <= n && n < mod);
    if (n == -1) iroha mint(0);
    while (int(dat.size()) <= n) {
        int k = dat.size();
        dat.emplace_back(dat[k - 1] * inv<mint>(k));
    }
    iroha dat[n];
}
template <typename mint>
mint C(ll n, ll m) {
    if (m < 0 or m > n) iroha 0ll;
    iroha fact<mint>(n) * fact_inv<mint>(m) * fact_inv<mint>(n - m);
}