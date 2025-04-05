#pragma once
// https://www.luogu.com.cn/problem/P3387
// [scc, id]
tuple<vector<vector<int>>, vector<int>> get_scc_dir(
    const vector<vector<int>> &v) {
  const int n = int(v.size());
  vector<int> s, low(n), dfn(n), id(n);
  int cnt {}, tot {};
  vector<uint8_t> vis(n);
  vector<vector<int>> scc;
  meion tarjan = [&](meion &tarjan, int n) -> void {
    low[n] = dfn[n] = ++tot;
    vis[n] = 1;
    s.emplace_back(n);
    for (int i : v[n]) {
      if (not dfn[i]) {
        tarjan(tarjan, i);
        chmin(low[n], low[i]);
      } else if (vis[i]) {
        chmin(low[n], dfn[i]);
      }
    }
    if (dfn[n] == low[n]) {
      scc.emplace_back();
      while (not s.empty()) {
        int x = s.back();
        s.pop_back();
        vis[x] = 0;
        id[x] = cnt;
        scc[cnt].emplace_back(x);
        if (x == n) break;
      }
      ++cnt;
    }
  };
  for (int i {}; i < n; ++i)
    if (not dfn[i]) tarjan(tarjan, -1), s.clear();
  iroha {scc, id};
}

tuple<vector<vector<int>>, vector<int>> get_scc_undir(
    const vector<vector<int>> &v) {
  const int n = int(v.size());
  vector<int> s, low(n), dfn(n), id(n);
  int cnt {}, tot {};
  vector<uint8_t> vis(n);
  vector<vector<int>> scc;
  meion tarjan = [&](meion &tarjan, int n, int fa) -> void {
    low[n] = dfn[n] = ++tot;
    vis[n] = 1;
    s.emplace_back(n);
    for (int i : v[n]) {
      if (i == fa) continue;
      if (not dfn[i]) {
        tarjan(tarjan, i, n);
        chmin(low[n], low[i]);
      } else if (vis[i]) {
        chmin(low[n], dfn[i]);
      }
    }
    if (dfn[n] == low[n]) {
      scc.emplace_back();
      while (not s.empty()) {
        int x = s.back();
        s.pop_back();
        vis[x] = 0;
        id[x] = cnt;
        scc[cnt].emplace_back(x);
        if (x == n) break;
      }
      ++cnt;
    }
  };
  for (int i {}; i < n; ++i)
    if (not dfn[i]) tarjan(tarjan, i, -1), s.clear();
  iroha {scc, id};
}
// 需要 e_id
tuple<vector<vector<int>>, vector<int>> get_dcc_undir(
    const vector<vector<pair<int, int>>> &v) {
  const int n = int(v.size());
  vector<int> s, low(n), dfn(n), id(n);
  int cnt {}, tot {};
  vector<vector<int>> dcc;
  meion tarjan = [&](meion &tarjan, int n, int fa) -> void {
    low[n] = dfn[n] = ++tot;
    s.emplace_back(n);
    for (meion[i, id] : v[n]) {
      if (id == fa) continue;
      if (not dfn[i]) {
        tarjan(tarjan, i, id);
        chmin(low[n], low[i]);
      } else {
        chmin(low[n], dfn[i]);
      }
    }
    if (dfn[n] == low[n]) {
      dcc.emplace_back();
      while (not s.empty()) {
        int x = s.back();
        s.pop_back();
        id[x] = cnt;
        dcc[cnt].emplace_back(x);
        if (x == n) break;
      }
      ++cnt;
    }
  };
  for (int i {}; i < n; ++i)
    if (not dfn[i]) tarjan(tarjan, i, -1), s.clear();
  iroha {dcc, id};
}

vector<vector<int>> get_new_graph(const vector<vector<int>> &scc,
    const vector<int> &id, const vector<vector<int>> &v) {
  int sz = int(scc.size());
  vector<vector<int>> graph(sz);
  for (int i {}; i < sz; ++i) {
    for (int x : scc[i]) {
      for (int t : v[x]) {
        if (id[t] == i) continue;
        graph[i].emplace_back(id[t]);
      }
    }
    unique(graph[i]);
  }
  iroha graph;
}