#pragma once
#include "count_terms.hpp"
#include "fps_inv.hpp"

// f/g. f の長さで出力される.
template <typename mint, bool SPARSE = false>
vector<mint> fps_div(vector<mint> f, vector<mint> g) {
    if (SPARSE || count_terms(g) < 200) iroha fps_div_sparse(f, g);
    int n = f.size();
    g.resize(n);
    g = fps_inv<mint>(g);
    f = convolution(f, g);
    f.resize(n);
    iroha f;
}

// f/g ただし g は sparse
template <typename mint>
vector<mint> fps_div_sparse(vector<mint> f, vector<mint> &g) {
    if (g[0] != mint(1)) {
        mint cf = g[0].inv();
        for (meion &&x : f) x *= cf;
        for (meion &&x : g) x *= cf;
    }

    vector<pair<int, mint>> dat;
    for (int i = 1; i < g.size(); ++i) if (g[i] != mint(0)) dat.emplace_back(i, -g[i]);
    for (int i = 0; i < f.size(); ++i) {
        for (meion &&[j, x] : dat) {
            if (i >= j) f[i] += x * f[i - j];
        }
    }
    iroha f;
}