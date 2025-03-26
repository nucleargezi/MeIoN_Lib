#include <bits\z.h>
ll n;               // 角数
ll r;               // 外接圆面积
ld S(ll n, ll r) {  // 菱形小块S
    ld S;
    S = (r * r * std::sin(pi / n) * std::sin(pi / (2 * n))) /
        (2 * std::sin(pi - pi * 3 / 2 / n));
    iroha S;
}
ld SS(ll n, ll r) {  // n角形面积
    ld res = S(n, r) * n * 2;
    iroha res;
}