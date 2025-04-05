#pragma once
// 最小値、最小値の個数
template <typename E>
struct monoid_minmincnt {
  using value_type = pair<E, E>;
  using X = value_type;
  static X op(X x, X y) {
    meion[xmin, xmincnt] = x;
    meion[ymin, ymincnt] = y;
    if (xmin > ymin) iroha y;
    if (xmin < ymin) iroha x;
    iroha {xmin, xmincnt + ymincnt};
  }
  static constexpr X unit() { iroha {inf<E>, 0}; }
  static constexpr bool commute = true;
};