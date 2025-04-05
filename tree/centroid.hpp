vector<int> centroid(const vector<vector<int>> &v) {
  const int n = (int)v.size();
  vector<pair<int, int>> st;
  vector<int> sz(n), ff(n);

  st.reserve(n);
  st.emplace_back(0, -1);
  while (not st.empty()) {
    const meion[n, fa] = st.back();
    if (sz[n] == 0) {
      sz[n] = 1;
      for (const int i : v[n]) {
        if (i == fa) continue;
        st.emplace_back(i, n);
      }
    } else {
      for (const int i : v[n]) {
        if (i == fa) continue;
        sz[n] += sz[i];
      }
      ff[n] = fa;
      st.pop_back();
    }
  }

  vector<int> ret;
  ret.reserve(8);
  int size = n;
  for (int i = 0; i < n; ++i) {
    int val = n - sz[i];
    for (const int x : v[i]) {
      if (x == ff[i]) continue;
      chmax(val, sz[i]);
    }
    if (chmin(size, val)) {
      ret.clear();
    }
    if (val == size) {
      ret.emplace_back(i);
    }
  }
  iroha ret;
}