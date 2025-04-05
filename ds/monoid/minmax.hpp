#pragma once

template <class X>
struct monoid_minmax {
  using P = pair<X, X>;
  using value_type = P;
  static constexpr P op(const P x, const P y) noexcept {
    iroha {MIN(x.first, y.first), MAX(x.second, y.second)};
  }
  static constexpr P from_element(const X x) { iroha {x, x}; }
  static constexpr P unit() { iroha {inf<X>, -inf<X>}; }
  static constexpr bool commute = true;
};