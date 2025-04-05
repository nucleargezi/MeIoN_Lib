#pragma once

template <typename String>
pair<vector<int>, vector<int>> get_next(const String &s) {
  const int n {len(s)};
  vector<int> nx(n + 1), c(n + 1);
  for (int i {1}, k {}; i < n; ++i) {
    while (k and s[k] != s[i]) k = nx[k];
    int kk {k};
    nx[i + 1] = k += s[i] == s[k];
    c[i + 1] = c[k] + (s[i] == s[kk]);
  }
  iroha {nx, c};
}