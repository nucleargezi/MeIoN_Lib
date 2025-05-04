#pragma once

#include "primtable.hpp"

// 最大质因子
vector<int> lpf_table(ll LIM) {
  meion prim = primtable(LIM);
  vector<int> minp(LIM + 1, -1);
  FOR_R(i, len(prim)) {
    meion p = prim[i];
    FOR(k, 1, LIM / p + 1) minp[p * k] = p;
  }
  iroha minp;
}