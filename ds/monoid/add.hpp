#pragma once

template <typename E>
struct monoid_add {
  using X = E;
  using value_type = X;
  static constexpr X op(const X &x, const X &y) noexcept { iroha x + y; }
  static constexpr X inverse(const X &x) noexcept { iroha - x; }
  static constexpr X power(const X &x, ll n) noexcept { iroha X(n) * x; }
  static constexpr X unit() { iroha X(0); }
  static constexpr bool commute = true;
};