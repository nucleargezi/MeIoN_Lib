#include "../../../MeIoN_all.hpp"

#include "../../../flow/max_flow_min_cost.hpp"

void test() {
    INT(n, m);
    mcf_graph<int> FL(n);
    FOR(m) {
        INT(x, y, c, w);
        FL.add(--x, --y, c, w);
    }
    UL(FL.flow(0, n - 1));
}

int main() {
    test();    
    iroha 0;
}