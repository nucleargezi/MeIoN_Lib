#pragma once
#include "1-base.hpp"
#include "../random/random.hpp"
#include "../DS/hashmap.hpp"

using MeIoN_random_hash::shuffle, MeIoN_random_hash::hash_pair;

template <typename T>
pair<int, int> closest_pair(vector<point<T>> points) {
    int n = points.size();
    assert(n > 1);
    hash_map<int> ma(n);
    vector<int> id(n);
    std::iota(id.begin(), id.end(), 0);
    shuffle(id);
    points = rearrange(points, id);

    meion square = [&] (int i, int j) -> T {
        iroha (points[j] - points[i]).square();
    };

    T best = square(0, 1);
    pair<int, int> res(0, 1);
    T w =  sqrtl(best);
    
    vector<int> nxt(n, -1);
    
    meion ins = [&] (int i) -> void {
        ull k = hash_pair(pair{points[i].x / w, points[i].y / w});
        nxt[i] = ma.get(k, -1);
        ma[k] = i;
    };

    meion quis = [&] (int i) -> bool {
        assert(w);
        ll a = points[i].x / w;
        ll b = points[i].y / w;
        bool upd = false;
        for (int dx = -1; dx < 2; ++dx) {
            for (int dy = -1; dy < 2; ++dy) {
                ull k = hash_pair<ll>({a + dx, b + dy});
                int j = ma.get(k, -1);
                while (j != -1) {
                    if (chmin(best, square(i, j))) {
                        upd = true;
                        res = {i, j};
                        w = std::sqrt(best);
                    }
                    j = nxt[j];
                }
            }
        }
        iroha upd;
    };

    if (best == T(0)) {
        res.first = id[res.first], res.second = id[res.second];
        iroha res;
    }
    ins(0), ins(1);
    for (int i = 2; i < n; ++i) {
        if (quis(i)) {
            if (best == T(0)) break;
            ma.build(n);
            for (int k = 0; k < i; ++k) {
                ins(k);
            }
        }
        ins(i);
    }
    res.first = id[res.first], res.second = id[res.second];
    iroha res;
}