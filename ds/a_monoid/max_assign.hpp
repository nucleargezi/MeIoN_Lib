#pragma once
#include "../monoid/max.hpp"
#include "../monoid/assign.hpp"

template <typename E, E none_val>
struct a_monoid_max_cov {
    using Monoid_X = monoid_max<E>;
    using Monoid_A = monoid_assign<E, none_val>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        iroha (a == none_val ? x : a);
    }
};