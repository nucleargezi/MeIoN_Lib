#pragma once

template <class X>
struct monoid_max {
  using value_type = X;
  static constexpr X op(const X &a, const X &b) noexcept {
    iroha std::max(a, b);
  }
  static constexpr X unit() { iroha std::numeric_limits<X>::lowest(); }
  static constexpr bool commute = true;
};