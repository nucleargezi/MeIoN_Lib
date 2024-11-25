#pragma once
#include "../monoid/add.hpp"
#include "../monoid/max.hpp"

template <typename E>
struct a_monoid_max_add {
    using Monoid_X = monoid_max<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (x == inf<E>) return x;
        return x + a;
    }
};