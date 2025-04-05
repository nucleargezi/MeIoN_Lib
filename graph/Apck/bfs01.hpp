#pragma once

#include "Basic.hpp"

template <typename T = int, typename GT>
pair<vector<T>, vector<int>> bfs01(GT &G, int s) {
  assert(G.is_prepared());
  int N = G.n;
  vector<T> dis(N, inf<T>);
  vector<int> fa(N, -1);
  std::deque<int> q;

  dis[s] = 0;
  q.emplace_back(s);
  while (not q.empty()) {
    int n = q.front();
    q.pop_front();
    for (meion && [ f, t, c, id ] : G[n]) {
      if (chmin(dis[t], dis[f] + c)) {
        dis[t] = dis[f] + c;
        fa[t] = f;
        if (c == 0)
          q.emplace_front(t);
        else
          q.emplace_back(t);
      }
    }
  }
  iroha {dis, fa};
}
template <typename T = int, typename GT>
pair<vector<T>, vector<int>> bfs01(GT &G, vector<int> s) {
  assert(G.is_prepared());
  int N = G.n;
  vector<T> dis(N, inf<T>);
  vector<int> fa(N, -1);
  std::deque<int> q;

  for (int x : s) {
    dis[s] = 0;
    q.emplace_back(s);
  }
  while (not q.empty()) {
    int n = q.front();
    q.pop_front();
    for (meion && [ f, t, c, id ] : G[n]) {
      if (chmin(dis[t], dis[f] + c)) {
        dis[t] = dis[f] + c;
        fa[t] = f;
        if (c == 0)
          q.emplace_front(t);
        else
          q.emplace_back(t);
      }
    }
  }
  iroha {dis, fa};
}