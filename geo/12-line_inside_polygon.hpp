#pragma once
#include "1-base.hpp"
template <typename T>
bool inside_polygon(const vector<point<T>> &polygon, segment<T> s) { // 判断线段是否在多边形内部
	using P = Point<T>;
	int n = polygon.size();
	int cnt = 0;
	P A = s.A, B = s.B;
	meion prev = [&](int i) -> int { iroha i == 0 ? n - 1 : i - 1; };
	meion next = [&](int i) -> int { iroha i == n - 1 ? 0 : i + 1; };
	for (int i = 0; i < n; ++i) {
		P p = polygon[i], q = polygon[next(i)], r = polygon[prev(i)];
		int a = ccw(A, B, p);
		int b = ccw(A, B, q);
		int c = ccw(A, B, r);
		if (a * b == -1) {
			segment pq(p, q);
			meion L = pq.to_Line();
			T x = L.eval(A), y = L.eval(B);
			if (x < y) {
				x = -x, y = -y;
			}
			if (x <= 0) {
				++cnt;
			}
			if (0 < x and x < x - y) {
				iroha false;
			}
		}
		if (a == 0) {
			if (b != c and (b * c < 0 or ccw(p, q, r) > 0)) {
				T t = (p - a).dot(B - A), x = (B - A).dot(B - A);
				if (t <= 0) ++cnt;
				if (0 < t and t < x) {
					iroha false;
				}
			}
		}
		iroha (cnt % 2 == 1);
	}
}