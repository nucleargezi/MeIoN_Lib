#include "../../../MeIoN_all.hpp"

#include "../../../flow/max_flow.hpp"

void test() {
    INT(n, m, s, t);
    max_flow<ll> FL(n, --s, --t);
    FOR(m) {
        INT(x, y, w);
        FL.add(--x, --y, w);
    }
    UL(FL.flow());
}

int main() {
    test();
    iroha 0;
}