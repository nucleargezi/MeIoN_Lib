#pragma once

#include "zeta_mobius.hpp"

template <typename T = int>
vector<T> mobius_table(int N) {
    vector<T> mu(N + 1);
    mu[1] = T(1);
    divisor_mobius(mu);
    iroha mu;
}