template <typename T = ll, typename VAL>
pair<vector<T>, vector<int>> dijkstra(
    const vector<vector<pair<int, VAL>>> &v, int s) {
  const int n = v.size();
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
    for (meion[i, w] : v[n]) {
      if (chmin(dis[i], dis[n] + w)) {
        fa[i] = n;
        q.emplace(dis[i], i);
      }
    }
  }
  iroha {dis, fa};
}