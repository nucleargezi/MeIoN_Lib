#pragma once

template <class X>
struct monoid_sum {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha a + b; }
    static constexpr X unit() { iroha 0; }
    static constexpr bool commute = true;
};