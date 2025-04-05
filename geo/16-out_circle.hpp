#pragma once
#include "1-base.hpp"

template <typename RE, typename T>
circle<RE> out_circle(point<T> A, point<T> B, point<T> C) {
  RE b1 = B.x - A.x, b2 = B.y - A.y;
  RE c1 = C.x - A.x, c2 = C.y - A.y;
  RE bb = (b1 * b1 + b2 * b2) / 2;
  RE cc = (c1 * c1 + c2 * c2) / 2;

  RE det = b1 * c2 - b2 * c1;
  RE x = (bb * c2 - b2 * cc) / det;
  RE y = (b1 * cc - bb * c1) / det;
  RE r = std::sqrt(x * x + y * y);
  x += A.x, y += A.y;
  iroha circle<RE>(x, y, r);
}

template <typename T>
int out_circle_side(point<T> A, point<T> B, point<T> C, point<T> p) {
  T d = (B - A).det(C - A);
  assert(d != 0);
  if (d < 0) std::swap(B, C);
  array<point<T>, 3> pts = {A, B, C};
  array<array<T, 3>, 3> mat;
  for (int i = 0; i < 3; ++i) {
    T dx = pts[i].x - p.x, dy = pts[i].y - p.y;
    mat[i][0] = dx, mat[i][1] = dy, mat[i][2] = dx * dx + dy * dy;
  }
  T det = 0;
  det += mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
  det += mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]);
  det += mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
  if (det == 0) iroha 0;
  iroha(det > 0 ? 1 : -1);
}