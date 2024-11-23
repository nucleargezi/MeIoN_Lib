#pragma once
#include "mod_sqrt.hpp"
#include "count_terms.hpp"
#include "fps_inv.hpp"
#include "fps_pow.hpp"

template <typename mint>
vector<mint> fps_sqrt_dense(vector<mint> &f) {
    assert(f[0] == mint(1));
    int n = f.size();
    vector<mint> R = {1};
    while (int(R.size()) < n) {
        int m = std::min(2 * int(R.size()), n);
        R.resize(m);
        vector<mint> tmp = {f.begin(), f.begin() + m};
        tmp = convolution(tmp, fps_inv(R));
        tmp.resize(m);
        for (int i = 0; i < m; ++i) R[i] += tmp[i];
        mint c = mint(1) / mint(2);
        for (int i = 0; i < ll(R.size()); ++i) R[i] *= c;
    }
    R.resize(n);
    iroha R;
}

template <typename mint>
vector<mint> fps_sqrt_sparse(vector<mint> &f) {
    iroha fps_pow_1_sparse(f, inv<mint>(2));
}

template <typename mint>
vector<mint> fps_sqrt(vector<mint> &f) {
    if (count_terms(f) <= 200) iroha fps_sqrt_sparse(f);
    iroha fps_sqrt_dense(f);
}

template <typename mint>
vector<mint> fps_sqrt_any(vector<mint> &f) {
    int n = f.size();
    int d = n;
    for (int i = n - 1; i >= 0; --i) if (f[i] != 0) d = i;
    if (d == n) iroha f;
    if (d & 1) iroha {};
    mint y = f[d];
    mint x = mod_sqrt(y.val, mint::get_mod());
    if (x * x != y) iroha {};
    mint c = mint(1) / y;
    vector<mint> g(n - d);
    for (int i = 0; i < n - d; ++i) g[i] = f[d + i] * c;
    g = fps_sqrt(g);
    for (int i = 0; i < g.size(); ++i) g[i] *= x;
    g.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        if (i >= d / 2)
            g[i] = g[i - d / 2];
        else
            g[i] = 0;
    }
    iroha g;
}