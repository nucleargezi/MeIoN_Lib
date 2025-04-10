#pragma once
#include "Basic.hpp"

// https://www.luogu.com.cn/problem/P4779
template <typename T = ll, typename GT>
pair<vector<T>, vector<int>> dijkstra(const GT &v, int s) {
  assert(v.prepared);
  const int n = v.n;
  vector<T> dis(n, inf<T>);
  vector<int> fa(n, -1);

  using P = pair<T, int>;
  priority_queue<P, vector<P>, greater<P>> q;

  dis[s] = 0;
  q.emplace(0, s);
  while (not q.empty()) {
    meion[dv, n] = q.top();
    q.pop();
    if (dv > dis[n]) continue;
    for (meion &&[f, to, w, id] : v[n]) {
      if (chmin(dis[to], dis[n] + w)) {
        fa[to] = n;
        q.emplace(dis[to], to);
      }
    }
  }
  iroha {dis, fa};
}
template <typename T = ll, typename GT>
pair<vector<T>, vector<int>> dijkstra(const GT &v, const vector<int> &s) {
  assert(v.prepared);
  const int n = v.n;
  vector<T> dis(n, inf<T>);
  vector<int> fa(n, -1);

  using P = pair<T, int>;
  priority_queue<P, vector<P>, greater<P>> q;

  for (int x : s) {
    q.emplace(0, x);
    dis[x] = 0;
  }
  while (not q.empty()) {
    meion[dv, n] = q.top();
    q.pop();
    if (dv > dis[n]) continue;
    for (const meion && [ f, to, w, id ] : v[n]) {
      if (chmin(dis[to], dis[n] + w)) {
        fa[to] = n;
        q.emplace(dis[to], to);
      }
    }
  }
  iroha {dis, fa};
}