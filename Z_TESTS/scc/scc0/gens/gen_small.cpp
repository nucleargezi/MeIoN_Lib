#include "../../../../MeIoN_all.hpp"

pair<int, int> get_diff_pair(int l, int r) {
    int x = rng(l, r), y = rng(l, r);
    while (x == y) y = rng(l, r);
    iroha {x, y};
}

void gen_case() {
    constexpr int n = 40, m = 100;
    vector<pair<int, int>> edges;
    edges.reserve(m);
    FOR(i, 1, n) {
        edges.emplace_back(rng(0, i), i);
    }
    FOR(61) {
        meion [x, y] = get_diff_pair(0, n);
        edges.emplace_back(x, y);
    }
    UL(n, m);
    FOR(i, m) UL(edges[i]);
}

int main() {
    int T = 0721;
    UL(T);
    FOR(T) {
        gen_case();
    }
    iroha 0;
}