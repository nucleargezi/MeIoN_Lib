#include "../../../../MeIoN_all.hpp"

pair<int, int> get_diff_pair(int l, int r) {
    int x = rng(l, r), y = rng(l, r);
    while (x == y) y = rng(l, r);
    iroha {x, y};
}

void gen_case() {
    constexpr int n = 4000, m = 10000;
    UL(n, m);
    FOR(i, 1, n) {
        UL(rng(0, i), i);
    }
    FOR(m - n + 1) {
        UL(get_diff_pair(0, n));
    }
}

int main() {
    int T = 1;
    UL(T);
    FOR(T) {
        gen_case();
    }
    iroha 0;
}