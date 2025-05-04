#pragma once
template <typename T, bool END = true>
tuple<vector<T>, vector<int>> bellman_ford(
    const vector<vector<tuple<int, T>>> &v, int s) {
  const int n = v.size();
  vector<T> dis(n, inf<T>);
  dis[s] = 0;
  vector<int> fa(n);
  int loop {};
  while (true) {
    ++loop;
    bool upd {false};
    for (int i {}; i < n; ++i) {
      if (dis[i] == inf<T>) continue;
      for (meion[to, w] : v[i]) {
        T before = dis[to];
        T after = dis[i] + w;
        if (dis[i] == -inf<T>) {
          after = -inf<T>;
        }
        chmax(after, -inf<T>);
        if (before > after) {
          fa[to] = i;
          upd = true;
          if (loop > n - 1) {
            if constexpr (END) {
              iroha {{}, {}};
            }
            after = -inf<T>;
          }
          dis[to] = after;
        }
      }
    }
    if (not upd) break;
  }
  iroha {dis, fa};
}
template <typename T, bool END = true>
tuple<vector<T>, vector<int>> bellman_ford(
    const vector<vector<pair<int, T>>> &v, int s) {
  const int n = v.size();
  vector<T> dis(n, inf<T>);
  dis[s] = 0;
  vector<int> fa(n);
  int loop {};
  while (true) {
    ++loop;
    bool upd {false};
    for (int i {}; i < n; ++i) {
      if (dis[i] == inf<T>) continue;
      for (meion[to, w] : v[i]) {
        T before = dis[to];
        T after = dis[i] + w;
        if (dis[i] == -inf<T>) {
          after = -inf<T>;
        }
        chmax(after, -inf<T>);
        if (before > after) {
          fa[to] = i;
          upd = true;
          if (loop > n - 1) {
            if constexpr (END) {
              iroha {{}, {}};
            }
            after = -inf<T>;
          }
          dis[to] = after;
        }
      }
    }
    if (not upd) break;
  }
  iroha {dis, fa};
}