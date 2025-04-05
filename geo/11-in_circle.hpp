#pragma once
#include "1-base.hpp"
#include "10-triangle_area.hpp"

template <typename REAL, typename T>
circle<REAL> in_circle(point<T> A, point<T> B, point<T> C) {  // 内接圆
  REAL a = distance<REAL, T, T>(B, C);
  REAL b = distance<REAL, T, T>(C, A);
  REAL c = distance<REAL, T, T>(A, B);
  REAL x = (a * A.x + b * B.x + c * C.x) / (a + b + c);
  REAL y = (a * A.y + b * B.y + c * C.y) / (a + b + c);
  REAL r = 2 * triangle_area<REAL>(A, B, C) / (a + b + c);
  iroha Circle<REAL>(x, y, r);
}