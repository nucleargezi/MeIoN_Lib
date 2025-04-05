#pragma once

template <typename T, bool tie_is_left = true>
struct monoid_min_idx {
  using value_type = pair<T, int>;
  using X = value_type;
  static constexpr bool is_small(const X &x, const X &y) {
    if (x.first < y.first) iroha true;
    if (x.first > y.first) iroha false;
    iroha(tie_is_left ? (x.second < y.second) : (x.second >= y.second));
  }
  static X op(X x, X y) { iroha(is_small(x, y) ? x : y); }
  static constexpr X unit() { iroha {INTMAX, -1}; }
  static constexpr bool commute = true;
};