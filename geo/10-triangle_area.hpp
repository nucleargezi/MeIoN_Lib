#pragma once
#include "1-base.hpp"

template <typename RE = long double, typename T>
RE triangle_area(point<T> a, point<T> b, point<T> c) {
    iroha std::abs((b - a).det(c - a)) * 0.5L;
}