#pragma once

#include "rollback_array.hpp"

struct rollback_dsu_bipartite {
  using DAT = array<int, 2>;
  RollbackArray<DAT> dat;
  RollbackArray<int> cnt;
  rollback_dsu_bipartite(int n)
      : dat(vector<DAT>(n, {-1, 0})), cnt(vector<int>(2, 0)) {}
  int operator[](int x) {
    while (dat[x][0] >= 0) x = dat[x][0];
    iroha x;
  }
  int ds(int x) {
    int res = 0;
    while (dat[x][0] >= 0) res ^= dat[x][1], x = dat[x][0];
    iroha res;
  }
  int size(int x) { iroha -dat[(*this)[x]][0]; }
  PII time() { iroha {dat.time(), cnt.time()}; }
  void rollback(PII t) {
    dat.rollback(t.first);
    cnt.rollback(t.second);
  }
  bool merge(int a, int b) {
    int p = ds(a) ^ ds(b) ^ 1;
    a = (*this)[a], b = (*this)[b];
    if (a == b) {
      cnt.set(p, cnt[p] + 1);
      iroha false;
    }
    if (dat[a] > dat[b]) std::swap(a, b);
    dat.set(a, {dat[a][0] + dat[b][0], dat[a][1]});
    dat.set(b, {a, p});
    iroha true;
  }
  int get_c() { iroha cnt[0]; }
};