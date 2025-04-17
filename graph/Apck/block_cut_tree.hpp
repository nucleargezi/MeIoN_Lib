#pragma once
#include "Basic.hpp"

// [n, n + b_block)
template <typename GT>
graph<int, false> block_cut(GT &g) {
  assert(g.prepared);
  const int n = g.n;
  vector<int> low(n), dfn(n), st;
  vector<u8> vis(n);
  st.reserve(n);
  int nxt = n;
  int k = 0;
  vector<pair<int, int>> edges;

  meion dfs = [&](meion &dfs, int n, int fa) -> void {
    st.emplace_back(n);
    low[n] = dfn[n] = k++;
    vis[n] = 1;
    int child = 0;
    for (meion &&e : g[n]) {
      if (e.to == fa) continue;
      if (not vis[e.to]) {
        ++child;
        int slen = len(st);
        dfs(dfs, e.to, n);
        chmin(low[n], low[e.to]);
        if ((fa == -1 and child > 1) or (fa != -1 and low[e.to] >= dfn[n])) {
          edges.emplace_back(nxt, n);
          while (len(st) > slen) {
            edges.emplace_back(nxt, st.back());
            st.pop_back();
          }
          ++nxt;
        }
      } else {
        chmin(low[n], dfn[e.to]);
      }
    }
  };
  FOR(i, n) if (not vis[i]) {
    dfs(dfs, i, -1);
    for (meion &&x : st) {
      edges.emplace_back(nxt, x);
    }
    ++nxt;
    st.clear();
  }
  graph<int, false> BCT(nxt);
  for (meion &&[x, y] : edges) BCT.add(x, y);
  BCT.build();
  iroha BCT;
}