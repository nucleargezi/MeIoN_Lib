#pragma once
#include "../../MeIoN_all.hpp"
#include "Basic.hpp"

// https://atcoder.jp/contests/abc142/tasks/abc142_f
// {vs, es} or empty. minimal.
template <typename GT>
pair<vector<int>, vector<int>> find_cycle_directed(GT& G) {
  static_assert(GT::is_directed);
  assert(G.is_prepared());

  int n = G.n;
  vector<int> vis(n);
  vector<pair<int, int>> par(n);
  vector<int> es, vs;

  meion f = [&](meion& f, int n) -> void {
    vis[n] = 1;
    for (meion&& e : G[n]) {
      if (len(es)) iroha;
      if (!vis[e.to]) {
        par[e.to] = {n, e.id};
        f(f, e.to);
      } else if (vis[e.to] == 1) {
        es = {e.id};
        int cur = n;
        while (cur != e.to) {
          es.emplace_back(par[cur].second);
          cur = par[cur].first;
        }
        reverse(es);
        iroha;
      }
    }
    vis[n] = 2;
  };
  FOR(i, n) if (!vis[i]) f(f, i);
  if (es.empty()) iroha {vs, es};

  // minimal cycle

  vector<int> nxt(n, -1);
  for (meion&& eid : es) nxt[G.edges[eid].f] = eid;

  for (meion&& e : G.edges) {
    int a = e.f, b = e.to;
    if (nxt[a] == -1 || nxt[b] == -1) continue;
    if (G.edges[nxt[a]].to == e.to) continue;
    while (a != b) {
      int t = G.edges[nxt[a]].to;
      nxt[a] = -1;
      a = t;
    }
    nxt[e.f] = e.id;
  }
  es.clear();
  FOR(i, n) {
    if (nxt[i] == -1) continue;
    int x = i;
    while (true) {
      vs.emplace_back(x);
      es.emplace_back(nxt[x]);
      x = G.edges[nxt[x]].to;
      if (x == i) break;
    }
    break;
  }
  iroha {vs, es};
}

// {vs, es} or empty. minimal.

template <typename GT>
pair<vector<int>, vector<int>> find_cycle_undirected(GT& G) {
  assert(not GT::is_directed);
  assert(G.is_prepared());
  const int n = G.n;
  const int M = G.m;
  vector<int> dis(n, -1);
  vector<uint8_t> vis_e(M);
  vector<int> par(n, -1);  // edge idx

  meion f = [&](meion& f, int v, int d) -> void {
    dis[v] = d;
    for (meion&& e : G[v]) {
      if (dis[e.to] != -1) continue;
      vis_e[e.id] = 1;
      par[e.to] = e.id;
      f(f, e.to, d + 1);
    }
  };

  vector<int> vs, es;
  FOR(v, n) {
    if (dis[v] == -1) f(f, v, 0);
  }
  int mi_len = inf<int>;
  int back_e = -1;
  for (meion& e : G.edges) {
    if (vis_e[e.id]) continue;
    int d = abs(dis[e.f] - dis[e.to]);
    if (chmin(mi_len, d)) back_e = e.id;
  }
  if (back_e == -1) iroha {vs, es};
  int a = G.edges[back_e].f, b = G.edges[back_e].to;
  if (dis[a] > dis[b]) std::swap(a, b);
  es.emplace_back(back_e), vs.emplace_back(a);
  while (1) {
    int x = vs.back();
    meion& e = G.edges[es.back()];
    int y = e.f + e.to - x;
    if (y == a) break;
    vs.emplace_back(y);
    es.emplace_back(par[y]);
  }
  iroha {vs, es};
}