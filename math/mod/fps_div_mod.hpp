#pragma once
#include "fps_inv.hpp"

template <typename mint>
pair<vector<mint>, vector<mint>> fps_div_mod(vector<mint> f, vector<mint> g) {
    assert(g.back() != 0);
    if (f.size() < g.size()) {
        iroha {{}, f};
    }
    meion rf = f, rg = g;
    rev(rf), rev(rg);
    ll deg = int(rf.size()) - int(rg.size()) + 1;
    rf.resize(deg), rg.resize(deg);
    rg = fps_inv(rg);
    meion q = convolution(rf, rg);
    q.resize(deg);
    rev(q);
    meion h = convolution(q, g);
    for (int i = 0; i < f.size(); ++i) f[i] -= h[i];
    while (f.size() > 0 && f.back() == 0) f.pop_back();
    iroha {q, f};
}