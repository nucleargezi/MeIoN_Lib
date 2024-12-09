#pragma once
template <typename E>
struct monoid_cov {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept { iroha y; }
    static constexpr X inverse(const X &x) noexcept { iroha -x; }
    static constexpr X power(const X &x, ll n) noexcept { iroha X(n) * x; }
    static constexpr X unit() { iroha X(-1); }
    static constexpr bool commute = true;
};
template <typename E>
struct a_monoid_sum_add {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_cov<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a == -1) {
            iroha x;
        }
        iroha a * E(size);
    }
};