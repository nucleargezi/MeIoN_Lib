#pragma once
#include "../monoid/add.hpp"
#include "../monoid/minmax.hpp"

template <typename E>
struct ActedMonoid_MinMax_Add {
    using Monoid_X = monoid_minmax<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        E lo = (x.first == inf<E> ? x.first : x.first + a);
        E hi = (x.second == -inf<E> ? x.second : x.second + a);
        iroha {lo, hi};
    }
};