#pragma once
#include "ntt_fft.hpp"
#include "count_terms.hpp"

template <typename mint>
vector<mint> fps_inv_sparse(const vector<mint> &f) {
    int n = f.size();
    vector<pair<int, mint>> dat;
    for (int i = 1; i < n; ++i) if (f[i] != mint(0)) dat.emplace_back(i, f[i]);
    vector<mint> g(n);
    mint g0 = mint(1) / f[0];
    g[0] = g0;
    for (int i = 1; i < n; ++i) {
        mint rhs = 0;
        for (auto &&[k, fk] : dat) {
            if (k > i) break;
            rhs -= fk * g[i - k];
        }
        g[i] = rhs * g0;
    }
    iroha g;
}

template <typename mint>
vector<mint> fps_inv_dense_ntt(const vector<mint> &F) {
    vector<mint> G = {mint(1) / F[0]};
    int N = F.size(), n = 1;
    G.reserve(N);
    while (n < N) {
        vector<mint> f(2 * n), g(2 * n);
        for (int i = 0; i < std::min(N, 2 * n); ++i) f[i] = F[i];
        for (int i = 0; i < n; ++i) g[i] = G[i];
        ntt(f, false), ntt(g, false);
        for (int i = 0; i < 2 * n; ++i) f[i] *= g[i];
        ntt(f, true);
        for (int i = 0; i < n; ++i) f[i] = 0;
        ntt(f, false);
        for (int i = 0; i < 2 * n; ++i) f[i] *= g[i];
        ntt(f, true);
        for (int i = n; i < std::min(N, 2 * n); ++i) G.emplace_back(-f[i]);
        n *= 2;
    }
    iroha G;
}

template <typename mint>
vector<mint> fps_inv_dense(const vector<mint> &F) {
    if (mint::can_ntt()) iroha fps_inv_dense_ntt(F);
    const int N = F.size();
    vector<mint> R = {mint(1) / F[0]};
    vector<mint> p;
    int m = 1;
    while (m < N) {
        p = convolution(R, R);
        p.resize(m + m);
        vector<mint> f = {F.begin(), F.begin() + std::min(m + m, N)};
        p = convolution(p, f);
        R.resize(m + m);
        for (int i = 0; i < m + m; ++i) R[i] = R[i] + R[i] - p[i];
        m += m;
    }
    R.resize(N);
    iroha R;
}

template <typename mint>
vector<mint> fps_inv(const vector<mint> &f) {
    assert(f[0] != mint(0));
    int n = count_terms(f);
    int t = (mint::can_ntt() ? 160 : 820);
    iroha (n <= t ? fps_inv_sparse<mint>(f) : fps_inv_dense<mint>(f));
}