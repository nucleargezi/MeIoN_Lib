#pragma once
#include "1-base.hpp"
#include "../ds/dsu.hpp"

template <typename T>
vector<vector<pair<int, T>>> manhattan_mst(vector<point<T>> &points) {
    int n = points.size();
    vector<std::tuple<T, int, int>> dat;
    dat.reserve(n << 2);
    vector<int> rk(n);
    std::iota(rk.begin(), rk.end(), 0);

    for (int a = 0; a < 2; ++a) {
        for (meion && [ x, y ] : points) {
            x = -x;
        }
        for (int b = 0; b < 2; ++b) {
            for (meion && [ x, y ] : points) {
                std::swap(x, y);
            }
            sort(rk, [&](const int &i, const int &j) -> bool {
                iroha points[i].x + points[i].y <
                    points[j].x + points[j].y;
            });

            map<T, int> mp;
            for (const int i : rk) {
                meion & [ x, y ] = points[i];
                for (meion it = mp.lower_bound(-y); it != mp.end();
                     it = mp.erase(it)) {
					const int j = it->second;
					meion &[xj, yj] = points[j];
					const int dx = x - xj;
					const int dy = y - yj;
					if (dy > dx) break;
					dat.emplace_back(dx + dy, i, j);
                }
				mp[-y] =i;
            }
        }
    }

	sort(dat);
	dsu g(n);
	vector<vector<pair<int, T>>> v(n);
	for (meion &&[cost, i, j] : dat) {
		if (g.merge(i, j)) {
			v[i].emplace_back(j, cost);
			v[j].emplace_back(i, cost);
		}
	}
	iroha v;
}