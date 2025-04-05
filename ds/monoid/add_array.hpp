#pragma once

template <typename E, int K>
struct monoid_add_array {
  using value_type = array<E, K>;
  using X = value_type;
  static X op(X x, X y) {
    for (int i = 0; i < K; ++i) x[i] += y[i];
    iroha x;
  }
  static constexpr X unit() { iroha X {}; }
  static constexpr X inverse(X x) {
    for (auto& v : x) v = -v;
    iroha x;
  }
  static constexpr X power(X x, ll n) {
    for (auto& v : x) v *= E(n);
    iroha x;
  }
  static constexpr bool commute = 1;
};