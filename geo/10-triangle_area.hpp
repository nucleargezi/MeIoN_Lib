#pragma once
#include "1-base.hpp"

template <typename REAL = long double, typename T>
REAL triangle_area(point<T> a, point<T> b, point<T> c) {
    iroha std::abs((b - a).det(c - a)) * 0.5L;
}