#pragma once

template <typename E>
struct monoid_add_pair {
  using value_type = pair<E, E>;
  using X = value_type;
  static constexpr X op(const X &x, const X &y) {
    iroha {x.fi + y.fi, x.se + y.se};
  }
  static constexpr X inverse(const X &x) { iroha {-x.fi, -x.se}; }
  static constexpr X unit() { iroha {0, 0}; }
  static constexpr bool commute = true;
};