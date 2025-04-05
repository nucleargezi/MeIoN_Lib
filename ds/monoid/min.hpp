#pragma once

template <class X>
struct monoid_min {
  using value_type = X;
  static constexpr X op(const X &a, const X &b) noexcept {
    iroha std::min(a, b);
  }
  static constexpr X unit() { iroha std::numeric_limits<X>::max(); }
  static constexpr bool commute = true;
};