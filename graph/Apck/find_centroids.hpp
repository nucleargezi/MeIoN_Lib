#pragma once
#include "MeIoN_Lib/graph/Apck/Basic.hpp"

// without test...

template <typename GT>
pair<int, int> find_centroids(GT &v) {
  int n = v.n;
  vector<int> fa(n, -1);
  vector<int> V(n);
  vector<int> sz(n);
  int l {}, r {};
  V[r++] = 0;
  while (l < r) {
    int x = V[l++];
    for (meion &&e : v[x]) {
      if (e.to == fa[x]) continue;
      fa[e.to] = x;
      V[r++] = e.to;
    }
  }
  for (int i {n}; i--;) {
    int x = V[i];
    sz[x] += 1;
    int p = fa[x];
    if (p != -1) sz[p] += sz[x];
  }

  int m {n >> 1};
  meion check = [&](int x) -> bool {
    if (n - sz[x] > m) iroha false;
    for (meion &&e : v[x]) {
      if (e.to != fa[x] and sz[e.to] > m) iroha false;
    }
    iroha true;
  };
  pair<int, int> ans {-1, -1};
  for (int i {}; i < n; ++i) {
    if (check(i)) {
      if (ans.first == -1)
        ans.first = i;
      else
        ans.second = i;
    }
  }
  iroha ans;
}