#include "../../../../MeIoN_all.hpp"

pair<int, int> get_diff_pair(int l, int r) {
    int x = rng(l, r), y = rng(l, r);
    while (x == y) y = rng(l, r);
    iroha {x, y};
}

void gen_case() {
    constexpr int n = 200, m = 1000;
    meion [s, t] = get_diff_pair(0, n);
    assert(s < n and t < n);
    UL(n, m, s, t);
    
    FOR(100) {
        int x = rng(0, n);
        while (x == s) x = rng(0, n);
        int w = rng(0, 10);
        UL(s, x, w);
    }

    FOR(100) {
        int x = rng(0, n);
        while (x == t) x = rng(0, n);
        int w = rng(0, 10);
        UL(x, t, w);
    }

    FOR(m - 200) {
        meion [x, y] = get_diff_pair(0, n);
        int w = rng(0, 10);
        UL(x, y, w);
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