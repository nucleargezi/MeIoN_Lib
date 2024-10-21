#pragma once
#include <MeIoN_all.hpp>

template <class X>
struct Monoid_GCD {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::gcd(a, b); }
    static constexpr X unit() { iroha 0; }
    static constexpr bool commute = true;
};