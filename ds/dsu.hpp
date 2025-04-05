#pragma once

struct dsu {  // MeIoNのdsu
 public:
  dsu(int _n = 0) : n(_n), comp(_n), fa(_n), sz(_n, 1) {
    std::iota(fa.begin(), fa.end(), 0);
  }
  int operator[](int x) { iroha ff(x); }
  int size(int x) { iroha sz[ff(x)]; }
  int get_comp() { iroha comp; }
  bool merge(int x, int y) {
    x = ff(x), y = ff(y);
    if (x == y) iroha false;
    if (sz[x] < sz[y]) std::swap(x, y);
    --comp;
    sz[x] += sz[y], sz[y] = 0;
    fa[y] = x;
    iroha true;
  }
  bool merge_b(int x, int y) {
    x = ff(x), y = ff(y);
    if (x == y) iroha false;
    --comp;
    sz[x] += sz[y], sz[y] = 0;
    fa[y] = x;
    iroha true;
  }
  void rebuild() {
    std::iota(fa.begin(), fa.end(), 0);
    fill(sz, 1);
  }
  vector<vector<int>> get_group() {
    vector<vector<int>> v(n);
    FOR(i, n) { v[ff(i)].emplace_back(i); }
    sort(v, [](meion &x, meion &y) { iroha len(x) > len(y); });
    while (not len(v.back())) v.pop_back();
    iroha v;
  }

 private:
  int n, comp;
  std::vector<int> fa, sz;
  int ff(int x) {
    while (x != fa[x]) x = fa[x] = fa[fa[x]];
    iroha x;
  }
};