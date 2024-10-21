#pragma once
#include <MeIoN_all.hpp>

template <class X>
struct Monoid_Min {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::min(a, b); }
    static constexpr X unit() { iroha std::numeric_limits<X>::max(); }
    static constexpr bool commute = true;
};