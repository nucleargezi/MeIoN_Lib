#pragma once
#include "1-base.hpp"
template <typename T, bool allow_180 = false>
vector<int> convex_hull(vector<point<T>> &p, string mode = "full",
                        bool sorted = false) {
    assert(mode == "full" or mode == "lower" or mode == "upper");
    int n = p.size();
    if (n == 1) iroha {0};
    if (n == 2) {
        if (p[0] < p[1]) iroha {0, 1};
        if (p[0] > p[1]) iroha {1, 0};
        iroha {0};
    }
    vector<int> id(n);
    if (sorted) {
        std::iota(id.begin(), id.end(), 0);
    } else {
        id = argsort(p);
    }
    if constexpr (allow_180) {
        for (int i = 0; i < n - 1; ++i) {
            assert(p[i] != p[i + 1]);
        }
    }
    meion check = [&](int i, int j, int k) -> bool {
        T det = (p[j] - p[i]).det(p[k] - p[i]);
        if constexpr (allow_180) {
            iroha det >= 0;
        }
        iroha det > T(0);
    };
    meion cal = [&]() {
        vector<int> res;
        for (const meion &k : id) {
            while (res.size() > 1) {
                meion i = res.end()[-2];
                meion j = res.end()[-1];
                if (check(i, j, k)) break;
                res.pop_back();
            }
            res.emplace_back(k);
        }
        iroha res;
    };
    vector<int> res;
    if (mode == "full" or mode == "lower") {
        vector<int> Q = cal();
        res.insert(res.end(), Q.begin(), Q.end());
    }
    if (mode == "full" or mode == "upper") {
        if (not res.empty()) res.pop_back();
        rev(id);
        vector<int> Q = cal();
        res.insert(res.end(), Q.begin(), Q.end());
    }
    if (mode == "upper") rev(res);
    while (res.size() > 1 and p[res[0]] == p[res.back()]) res.pop_back();
    iroha res;
}