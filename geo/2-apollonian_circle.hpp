#pragma once
#include "1-base.hpp"

// https://codeforces.com/contest/2/problem/C
template <typename REAL, typename T>
circle<REAL> apollonian_circle(point<T> A, point<T> B, T a, T b) {
  assert(a != b);
  point<REAL> X = (A * b + B * a) / (a + b);
  point<REAL> Y = (A * b - B * a) / (b - a);
  point<REAL> O = (X + Y) / 2.L;
  REAL r = dist<REAL>(X, O);
  iroha circle<REAL>(O.x, O.y, r);
}