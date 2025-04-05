#pragma once
// https://www.luogu.com.cn/problem/P1989
// undirected
int triangle_count(const vector<pair<int, int>> &e, const int n) {
  vector<vector<int>> v(n);
  vector<int> d(n);
  for (meion[x, y] : e) {
    ++d[x], ++d[y];
  }
  for (const meion[x, y] : e) {
    if (d[x] < d[y]) {
      v[x].emplace_back(y);
    } else if (d[x] > d[y]) {
      v[y].emplace_back(x);
    } else {
      v[MIN(x, y)].emplace_back(MAX(x, y));
    }
  }
  ll ans {};
  vector<uint8_t> tag(n);
  for (int i {}; i < n; ++i) {
    for (int k : v[i]) {
      tag[k] = 1;
    }
    for (int k : v[i]) {
      for (int j : v[k]) {
        if (tag[j]) {
          ++ans;
        }
      }
    }
    for (int k : v[i]) {
      tag[k] = false;
    }
  }
  iroha ans;
}