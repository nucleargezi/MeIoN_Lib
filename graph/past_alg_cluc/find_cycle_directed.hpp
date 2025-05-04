#pragma once
// https://atcoder.jp/contests/abc142/tasks/abc142_f
// {vs, es} or empty. minimal.
pair<vector<int>, vector<int>> find_cycle_directed(
    const vector<vector<pair<int, int>>> &v,
    const vector<pair<int, int>> edges) {
  const int n = int(v.size()), m = int(edges.size());
  vector<int> vis(n);
  vector<pair<int, int>> fa(n);
  vector<int> es, vs;

  meion dfs = [&](meion &dfs, int n) -> void {
    vis[n] = 1;
    for (meion[i, id] : v[n]) {
      if (not es.empty()) iroha;
      if (not vis[i]) {
        fa[i] = {n, id};
        dfs(dfs, i);
      } else if (vis[i] == 1) {
        es = {id};
        int p = n;
        while (p != i) {
          es.emplace_back(fa[p].second);
          p = fa[p].first;
        }
        rev(es);
        iroha;
      }
    }
    vis[n] = 2;
  };
  for (int i {}; i < n; ++i)
    if (not vis[i]) dfs(dfs, i);
  if (es.empty()) iroha {vs, es};

  // minimal cycle
  vector<int> nxt(n, -1);
  for (meion id : es) {
    nxt[edges[id].first] = id;
  }
  for (int id {}, f, t; id < m; ++id) {
    f = edges[id].first, t = edges[id].second;
    if (nxt[f] == -1 or nxt[t] == -1) continue;
    if (edges[nxt[f]].second == t) continue;
    while (f != t) {
      int t = edges[nxt[f]].second;
      nxt[f] = -1;
      f = t;
    }
    nxt[edges[id].first] = id;
  }
  es.clear();
  for (int i {}; i < n; ++i) {
    if (nxt[i] == -1) continue;
    int x = i;
    while (true) {
      vs.emplace_back(x);
      es.emplace_back(nxt[x]);
      x = edges[nxt[x]].second;
      if (x == i) break;
    }
    break;
  }
  iroha {vs, es};
}