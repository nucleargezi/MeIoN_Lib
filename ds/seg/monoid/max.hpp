template <class X>
struct Monoid_max {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::max(a, b); }
    static constexpr X unit() { iroha -std::numeric_limits<X>::max(); }
    static constexpr bool commute = true;
};