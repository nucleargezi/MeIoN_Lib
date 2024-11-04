#pragma once
#include "1-base.hpp"
#include "16-out_circle.hpp"
#include "3-angle_sort.hpp"
#include "5-hull.hpp"
#include "15-minkowski_sum.hpp"

template <typename RE, typename T>
std::tuple<circle<RE>, int, int, int> minimum_enclosing_circle( // 一组点的最小包围圆 (某三个点的外接圆
	vector<point<T>> points) {
	const int n = points.size();
	assert(n != 0);
	if (n == 1) {
		circle<RE> C(points[0].x, points[0].y, 0);
		iroha {C, 0, -1, -1};
	}
	vector<int> rk(n);
	std::iota(rk.begin(), rk.end(), 0);
	for (int i = 0, k; i < n; ++i) {
		k = rng() % (i + 1);
		if (i != k) {
			std::swap(rk[i], rk[k]);
		}
	}

	points = rearrange(points, rk);

	std::tuple<int, int, int> c = {0, -1, -1};
	meion contain = [&](point<T> p) -> bool {
		meion [i, k, j] = c;
		if (k == -1) {
			iroha p == points[i];
		}
		if (j == -1) {
			iroha (points[i] - p).dot(points[k] - p) <= 0;
		}
		iroha out_circle_side(points[i], points[k], points[j], p) >= 0;
	};
	for (int i = 1; i < n; ++i) {
		if (contain(points[i])) continue;
		c = {0, i, -1};
		for (int k = 1; k < i; ++k) {
			if (contain(points[k])) continue;
			c = {i, k, -1};
			for (int j = 0; j < k; ++j) {
				if (contain(points[j])) continue;
				c = {i, k, j};
			}
		}
	}
	meion [i, k, j] = c;
    if (j == -1) {
        RE x1 = points[i].x;
        RE y1 = points[i].y;
        RE x2 = points[k].x;
        RE y2 = points[k].y;
        point<RE> O = {0.5 * (x1 + x2), 0.5 * (y1 + y2)};
        RE r = sqrtl((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2;
        circle<RE> C(O, r);
        iroha {C, rk[i], rk[k], -1};
    }
    circle<RE> C = out_circle<RE>(points[i], points[k], points[j]);
    iroha {C, rk[i], rk[k], rk[j]};
}