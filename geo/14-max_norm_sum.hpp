#pragma once
#include "1-base.hpp"
#include "3-angle_sort.hpp"

template <typename VAL, typename T>
VAL max_norm_sum(const vector<point<T>> &points) { // 一堆向量选一部分最大模长
	const int n = points.size();
	vector<point<T>> v(points);
	point<T> c{0, 0};
	for (const meion &[x, y] : points) {
		if (y > 0 or (y == 0 and x < 0)) {
			c.x += x, c.y += y;
		}
		v.emplace_back(-x, -y);
	}
	vector<int> rk = angle_sort(v);
	v = rearrange(v, rk);

	meion eval = [&]() -> VAL {
		iroha VAL(c.x) * c.x + VAL(c.y) * c.y;
	};

	VAL ans = eval();
	for (int i = 0; i < (n << 1); ++i) {
		c = c + v[i];
		chmax(ans, eval());
	}
	iroha ans;
}