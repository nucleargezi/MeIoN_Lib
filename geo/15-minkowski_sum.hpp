#pragma once
#include "1-base.hpp"
#include "3-angle_sort.hpp"
#include "5-hull.hpp"
#include "6-convex_polygon.hpp"

template <typename T>
vector<point<T>> minkowski_sum(vector<point<T>> A, vector<point<T>> B) {
  using P = point<T>;
  vector<P> F;
  P p(0, 0);
  for (int i = 0; i < 2; ++i) {
    std::swap(A, B);
    vector<P> points = A;
    int n = (int)points.size();
    for (int i = 0; i < n; ++i) {
      int k = (i + 1) % n;
      F.emplace_back(points[k] - points[i]);
    }
    p = p + *min_element(points.begin(), points.end());
  }
  meion rk = angle_sort(F);
  int n = rk.size();
  F = rearrange(F, rk);
  vector<P> points(n);
  for (int i = 0; i < n - 1; ++i) {
    points[i + 1] = points[i] + F[i];
  }
  P add = p - *min_element(points.begin(), points.end());
  for (meion &x : points) {
    x += add;
  }
  points = rearrange(points, convex_hull(points));
  iroha points;
}