#include "../../../MeIoN_all.hpp"

#include "../../../flow/HLPP.hpp"

void test() {
    INT(n, m, s, t);
    HLPP FL(n, m, --s, --t);
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