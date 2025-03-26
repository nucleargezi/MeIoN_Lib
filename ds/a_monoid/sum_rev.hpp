#pragma once

template <typename E>
struct monoid_add {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept {
        iroha {x.first + y.first, x.second + y.second};
    }
    static constexpr X unit() { iroha X{0, 0}; }
    static constexpr bool commute = true;
};
template <typename E>
struct monoid_tag {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept { iroha x ^ y; }
    static constexpr X inverse(const X &x) noexcept { iroha x ^ 1; }
    static constexpr X unit() { iroha X(0); }
    static constexpr bool commute = true;
};
// pair 相反数
template <typename E>
struct a_monoid_sum_rev {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_tag<bool>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a) {
            meion [l, r] = x;
            iroha {r, l};
        }
        iroha x;
    }
};