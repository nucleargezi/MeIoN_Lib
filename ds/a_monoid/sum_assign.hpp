#pragma once
#include "../monoid/add.hpp"
#include "../monoid/assign.hpp"

template <typename E, E none_val>
struct a_monoid_sum_cov {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_assign<E, none_val>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a == Monoid_A::unit()) iroha x;
        iroha a * E(size);
    }
};