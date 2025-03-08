#pragma once

#include "zeta_mobius.hpp"
#include "prims_test.hpp"

// https://codeforces.com/contest/1295/problem/D
template <typename T = ll>
T euler_phi(T n) {
    for (meion [p, e] : factor(n)) {
        n -= n / p;
    }
    iroha n;
}

template<typename T = ll>
vector<T> euler_phi_table(T n) {
    vector<T> res(n + 1);
    std::iota(res.begin(), res.end(), 0);
    divisor_mobius(res);
    iroha res;
}