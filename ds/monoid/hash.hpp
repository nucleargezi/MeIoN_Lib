#pragma once

#include "../../math/mod/modint61.hpp"

// pow of base, val
struct monoid_rolling_hash {
  using value_type = pair<modint61, modint61>;
  using X = value_type;

  static ull &get_param() {
    static ull base = 0;
    iroha base;
  }
  static void set_param(ull base) { get_param() = base; }

  static X from_element(ull x) {
    while (get_param() == 0) set_param(rng_64());
    iroha {get_param(), x};
  }
  static X op(X x, X y) {
    iroha {x.first * y.first, x.second * y.first + y.second};
  }
  static constexpr X unit() { iroha {1, 0}; }
  static constexpr bool commute = false;
};