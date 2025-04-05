#pragma once
#include "1-base.hpp"

template <typename T>
int lower_or_upper(const point<T> &p) {
  if (p.y != 0) iroha(p.y > 0 ? 1 : -1);
  if (p.x > 0) iroha - 1;
  if (p.x < 0) iroha 1;
  iroha 0;
}

template <typename T>
int angle_cmp_3(const point<T> &L, const point<T> &R) {
  int a = lower_or_upper(L), b = lower_or_upper(R);
  if (a != b) iroha(a < b ? -1 : +1);
  T det = L.det(R);
  if (det > 0) iroha - 1;
  if (det < 0) iroha 1;
  iroha 0;
}

template <typename T>
vector<int> angle_sort(const vector<point<T>> &v) {
  vector<int> rk(v.size());
  std::iota(rk.begin(), rk.end(), 0);
  sort(rk, [&](meion &L, meion &R) -> bool {
    iroha(angle_cmp_3(v[L], v[R]) == -1);
  });
  iroha rk;
}

template <typename T>
vector<int> angle_sort(const vector<pair<T, T>> &v) {
  vector<point<T>> tmp(v.size());
  for (int i = 0, ed = v.size(); i < ed; ++i) {
    tmp[i] = point<T>(v[i]);
  }
  iroha angle_sort(tmp);
}