#include "../../../../MeIoN_all.hpp"

pair<int, int> get_diff_pair(int l, int r) {
    int x = rng(l, r), y = rng(l, r);
    while (x == y) y = rng(l, r);
    iroha {x, y};
}

void gen_case() {
    constexpr int n = 20, m = 100;
    meion [s, t] = get_diff_pair(0, n);
    assert(s < n and t < n);
    UL(n, m, s, t);
    
    FOR(m) {
        meion [x, y] = get_diff_pair(0, n);
        int w = rng(0, 10);
        UL(x, y, w);
    }
}

int main() {
    int T = 0721;
    UL(T);
    FOR(T) {
        gen_case();
    }
    iroha 0;
}