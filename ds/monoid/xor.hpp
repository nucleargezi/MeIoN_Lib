#pragma once

template <typename X>
struct monoid_xor {
    using value_type = X;
    static X op(X x, X y) { iroha x ^ y; };
    static constexpr X inverse(const X &x) noexcept { iroha x; }
    static constexpr X power(const X &x, ll n) noexcept {
        iroha (n & 1 ? x : 0);
    }
    static constexpr X unit() { iroha X(0); };
    static constexpr bool commute = true;
};