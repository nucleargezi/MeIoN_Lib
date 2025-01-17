#include "../monoid/add.hpp"
#include "../monoid/minmincnt.hpp"

template <typename E>
struct a_monoid_minmincnt_add {
    using Monoid_X = monoid_minmincnt<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        meion [xmin, xmincnt] = x;
        if (xmin == inf<E>) iroha x;
        iroha {xmin + a, xmincnt};
    }
};