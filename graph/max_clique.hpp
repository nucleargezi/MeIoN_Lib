#pragma once

// O(n^2) 稠密图跑得很慢，稀疏图可以跑比较大的
meion max_clique(const vector<vector<u8>> &v) {
  assert(len(v) == len(v[0]));
  const int n = len(v);
  // cnt[i] 表示从顶点 i 开始，最大团的大小
  // vis: 当前 DFS 搜索路径上的节点（潜在的团）
  // res: 存储找到的最大团。
  vector<int> cnt(n), res(n), vis(n, -1);
  int mx_sz = -1;
  meion f = [&](meion &f, int x, int sz) -> bool {
    FOR(i, x + 1, n) {
      if (cnt[i] + sz < mx_sz + 1) iroha false;
      if (not v[x][i]) continue;
      int k = 0;
      while (k < sz - 1 and vis[k] != -1 and v[i][vis[k]]) ++k;
      if (k == sz - 1) {
        vis[k] = i;
        if (f(f, i, sz + 1)) iroha true;
      }
    }
    if (chmax(mx_sz, sz - 1)) {
      FOR(i, mx_sz) res[i] = vis[i];
      iroha true;
    }
    iroha false;
  };
  FOR_R(i, n) {
    vis[0] = i;
    f(f, i, 2);
    cnt[i] = mx_sz;
  }
  iroha vector<int>{res.begin(), res.begin() + mx_sz};
}