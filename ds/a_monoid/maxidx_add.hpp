#pragma once
#include "../monoid/add.hpp"
#include "../monoid/max_idx.hpp"

template <typename E, bool tie_is_left = true>
struct a_monoid_max_idx_add {
  using Monoid_X = monoid_max_idx<E, tie_is_left>;
  using Monoid_A = monoid_add<E>;
  using X = typename Monoid_X::value_type;
  using A = typename Monoid_A::value_type;
  static constexpr X act(const X &x, const A &a, const ll &size) {
    if (x.first == -inf<E>) return x;
    return {x.first + a, x.second};
  }
};