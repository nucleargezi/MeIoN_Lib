#pragma once
#include "../monoid/minmax.hpp"

template <typename E>
struct monoid_tag {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept { iroha x ^ y; }
    static constexpr X inverse(const X &x) noexcept { iroha x ^ 1; }
    static constexpr X unit() { iroha X(0); }
    static constexpr bool commute = true;
};
// 相反数
template <typename E>
struct a_monoid_minmax_rev {
    using Monoid_X = monoid_minmax<E>;
    using Monoid_A = monoid_tag<bool>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a) {
            meion [min, max] = x;
            iroha {-max, -min};
        }
        iroha x;
    }
};