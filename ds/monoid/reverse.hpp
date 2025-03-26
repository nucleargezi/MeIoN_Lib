#pragma once

template <class monoid>
struct monoid_reverse {
    using value_type = typename monoid::value_type;
    using X = value_type;
    static constexpr X op(const X &x, const X &y) { return monoid::op(y, x); }
    static constexpr X unit() { return monoid::unit(); }
    static const bool commute = monoid::commute;
};