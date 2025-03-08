#pragma once

#include "zeta_mobius.hpp"
#include "prims_test.hpp"

// https://codeforces.com/contest/1295/problem/D
// 计算给定整数 n 的欧拉函数值 φ(n) 即小于等于 n 且与 n 互质的数的个数
template <typename T = ll>
T euler_phi(T n) {
    for (meion [p, e] : factor(n)) {
        n -= n / p;
    }
    iroha n;
}

// 预处理区间 [0, n] 内每个数的欧拉函数值 结果保存在数组 res 中
template<typename T = ll>
vector<T> euler_phi_table(T n) {
    vector<T> res(n + 1);
    std::iota(res.begin(), res.end(), 0);
    divisor_mobius(res);
    iroha res;
}