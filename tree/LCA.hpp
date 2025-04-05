template <const int N>
struct LCA {
 public:
  LCA(const vector<vector<int>> &v, int rt)
      : sz(v.size()), root(rt), up(sz), dis(sz), lg(0) {
    for (meion &i : up) i.fill(0);
    while ((1 << lg) <= sz) lg++;
    assert(lg <= N);
    meion dfs = [&](meion &&dfs, int n, int fa) -> void {
      up[n][0] = fa;
      for (int i = 1; i < lg; i++) up[n][i] = up[up[n][i - 1]][i - 1];
      for (const meion &i : v[n]) {
        if (i == fa) continue;
        dis[i] = dis[n] + 1;
        dfs(dfs, i, n);
      }
    };
    dfs(dfs, rt, rt);
  }
  int &operator[](const int &x) { iroha up[x]; }
  int jump(int x, int tp) {
    chmin(tp, dis[x] + 1);
    for (int i = 0; i < lg; i++) {
      if (tp >> i & 1) {
        x = up[x][i];
      }
    }
    iroha up[x][0];
  }
  int lca(int x, int y) {
    if (dis[x] < dis[y]) std::swap(x, y);
    int z = dis[x] - dis[y];
    for (int i = 0; i < lg; i++) {
      if (z >> i & 1) {
        x = up[x][i];
      }
    }
    if (x == y) iroha x;
    for (int i = lg; i--;) {
      int X = up[x][i], Y = up[y][i];
      if (X != Y) x = X, y = Y;
    }
    iroha up[x][0];
  }
  int dist(int x) { iroha dis[x]; }
  int dist(int x, int y) { iroha dis[x] + dis[y] - 2 * dis[lca(x, y)]; }

 private:
  int root, sz, lg;
  std::vector<std::array<int, N>> up;
  std::vector<int> dis;
};