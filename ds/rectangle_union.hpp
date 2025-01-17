#pragma once
#include "seg/lazy_seg_base.hpp"
#include "a_monoid/minmincnt_add.hpp"

template <typename XY = int>
struct rectangle_union {
    using RECT = tuple<XY, XY, XY, XY>;
    vector<RECT> rectangles;
    vector<XY> X, Y;

    void add(XY xl, XY yl, XY xr, XY yr) {
        assert(xl < xr && yl < yr);
        X.emplace_back(xl), X.emplace_back(xr), Y.emplace_back(yl),
            Y.emplace_back(yr);
        rectangles.emplace_back(xl, xr, yl, yr);
    }

    template <typename ANS_TYPE = ll>
    ANS_TYPE calc() {
        int N = X.size();
        vector<int> ord_x = argsort(X);
        vector<int> ord_y = argsort(Y);
        vector<int> rk_y(N);
        for (int i{}; i < N; ++i) rk_y[ord_y[i]] = i;
        X = rearrange(X, ord_x);
        Y = rearrange(Y, ord_y);

        using AM = a_monoid_minmincnt_add<XY>;
        lazy_seg<AM> seg(
            N - 1, [&](int i) -> pair<XY, XY> { iroha {0, Y[i + 1] - Y[i]}; });

        ANS_TYPE ANS = 0;
        XY total = Y.back() - Y[0];
        for (int i{}; i < N - 1; ++i) {
            int k = ord_x[i] / 2;
            int a = (ord_x[i] & 1 ? -1 : 1);
            seg.apply(rk_y[2 * k], rk_y[2 * k + 1], a);
            auto [min, mincnt] = seg.prod_all();
            ANS_TYPE dy = total - (min == 0 ? mincnt : 0);
            ANS_TYPE dx = X[i + 1] - X[i];
            ANS += dx * dy;
        }
        iroha ANS;
    }
};