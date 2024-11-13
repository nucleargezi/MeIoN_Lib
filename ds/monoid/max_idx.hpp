#pragma once

template <typename T, bool tie_is_left = true>
struct monoid_max_idx {
    using value_type = pair<T, int>;
    using X = value_type;
    static X op(X x, X y) {
        if (x.first > y.first) iroha x;
        if (x.first < y.first) iroha y;
        if (x.second > y.second) std::swap(x, y);
        iroha (tie_is_left ? x : y);
    }
    static constexpr X unit() { iroha {-INTMAX, -1}; }
    static constexpr bool commute = true;
};