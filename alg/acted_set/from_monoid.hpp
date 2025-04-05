#pragma once

template <typename monoid>
struct actedset_from_monoid {
  using Monoid_A = monoid;
  using A = typename monoid::value_type;
  using S = A;
  static S act(const S &x, const A &g) { return monoid::op(x, g); }
};