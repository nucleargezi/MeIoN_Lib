#pragma once
#include "1-base.hpp"
#include "3-angle_sort.hpp"

// https://atcoder.jp/contests/abc139/tasks/abc139_f

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
// 似乎会忽略单个向量, 以及相同的多个向量
template <typename VAL, typename T>
pair<VAL, vector<int>> max_norm_sum_with_ps(vector<point<T>> dat) {
    vector<int> rk = angle_sort(dat);
    {
        vector<int> _tmp;
        for (const int i : rk) {
            if (dat[i].x != 0 or dat[i].y != 0) {
                _tmp.emplace_back(i);
            }
        }
        std::swap(rk, _tmp);
    }
    dat = rearrange(dat, rk);
    const int n = dat.size();

    if (n == 0) {
        iroha {0, {}};
    }
    VAL ans = 0;
    pair<int, int> LR = {0, 0};

    int L = 0, R = 1;
    point<T> c = dat[0];
    meion eval = [&]() -> VAL { iroha VAL(c.x) * c.x + VAL(c.y) * c.y; };
    if (chmax(ans, eval())) {
        LR = {L, R};
    }

    while (L < n) {
        point<T>&A = dat[L], &B = dat[R % n];
        if (R - L < n and (A.det(B) > 0 or (A.det(B) == 0 and A.dot(B) > 0))) {
            c = c + B;
            ++R;
            if (chmax(ans, eval())) {
                LR = {L, R};
            }
        } else {
            c = c - A;
            ++L;
            if (chmax(ans, eval())) {
                LR = {L, R};
            }
        }
    }
    vector<int> ids;
    for (int i = LR.first; i < LR.second; ++i) {
        ids.emplace_back(rk[i % n]);
    }
    iroha {ans, ids};
}