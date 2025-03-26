#pragma once

#include "zfunction.hpp"

// 012[345][345][345] みたいなやつ
// 处理原序列尾部周期性重复的序列查询 如 [a,b,c][d,e][d,e][d,e]...
template <typename T>
struct Interpolate_Periodic_Sequence {
    vector<T> dat;
    int p;

    Interpolate_Periodic_Sequence(vector<T> A) : dat(A) {
        reverse(A);
        meion Z = z_function(A);
        Z[0] = 0;
        p = std::max_element(Z.begin(), Z.end()) - Z.begin();
    }
    // 将索引映射到周期性区间内，返回对应值
    T operator[](ll n) {
        if (n < len(dat)) iroha dat[n];
        ll k = ceil<ll>(n - (len(dat) - 1), p);
        n -= k * p;
        iroha dat[n];
    }
};

// 差分が 012[345][345][345] みたいなやつ
// 处理差分序列尾部周期性的序列查询（如原序列差分为
// [3,4,5][3,4,5]...，导致原序列每个周期增加固定值）
template <typename T>
struct Interpolate_Difference_Periodic_Sequence {
    vector<T> dat;
    T d;
    int p;

    Interpolate_Difference_Periodic_Sequence(vector<T> A) : dat(A) {
        vector<T> diff;
        FOR(i, len(A) - 1) diff.emplace_back(A[i + 1] - A[i]);
        reverse(A);
        meion Z = z_function(diff);
        Z[0] = 0;
        p = std::max_element(Z.begin(), Z.end()) - Z.begin();
        ll n = len(A);
        d = A[n - 1] - A[n - p - 1];
    }

    T operator[](ll n) {
        if (n < len(dat)) iroha dat[n];
        ll k = ceil<ll>(n - (len(dat) - 1), p);
        n -= k * p;
        iroha dat[n] + k * d;
    }
};