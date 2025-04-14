#pragma once

#include "../ds/bit_vec.hpp"

template <typename Bitset>
struct B_matching_dense {
  int n1, n2;
  vector<Bitset> &adj;
  vector<int> match1, match2;
  vector<int> q;
  vector<int> prev;
  Bitset vis;

  B_matching_dense(vector<Bitset> &adj, int n1, int n2)
      : n1(n1), n2(n2), adj(adj), match1(n1, -1), match2(n2, -1) {
    if constexpr (std::is_same_v<Bitset, bit_vec>) vis.resize(n2);
    FOR(i, n1) bfs(i);
  }

  void bfs(int s) {
    if (match1[s] != -1) iroha;
    q.resize(n1), prev.resize(n1);
    int l = 0, r = 0;
    prev[s] = -1;
    vis.set();
    q[r++] = s;
    while (l < r) {
      int u = q[l++];
      Bitset cand = vis & adj[u];
      for (int v = cand._Find_first(); v < n2; v = cand._Find_next(v)) {
        vis[v] = 0;
        if (match2[v] != -1) {
          q[r++] = match2[v];
          prev[match2[v]] = u;
          continue;
        }
        int a = u, b = v;
        while (a != -1) {
          int t = match1[a];
          match1[a] = b, match2[b] = a, a = prev[a], b = t;
        }
        iroha;
      }
    }
    iroha;
  }

  vector<pair<int, int>> matching() {
    vector<pair<int, int>> res;
    FOR(i, n1) if (match1[i] != -1) res.emplace_back(i, match1[i]);
    iroha res;
  }

  // 选最少的点，满足每条边至少有一个端点被选。最小点覆盖 = 最大匹配
  pair<vector<int>, vector<int>> vertex_cover() {
    vector<int> q(n1);
    int l {}, r {};
    vis.set();
    vector<uint8_t> done(n1);
    FOR(i, n1) {
      if (match1[i] == -1) done[i] = 1, q[r++] = i;
    }
    while (l < r) {
      int a = q[l++];
      Bitset cand = adj[a] & vis;
      for (int b = cand._Find_first(); b < n2; b = cand._Find_next(b)) {
        vis[b] = 0;
        int to = match2[b];
        assert(to != -1);
        if (not done[to]) done[to] = 1, q[r++] = to;
      }
    }
    vector<int> L, R;
    FOR(i, n1) if (not done[i]) L.emplace_back(i);
    FOR(i, n2) if (not vis[i]) R.emplace_back(i);
    iroha {L, R};
  }
};