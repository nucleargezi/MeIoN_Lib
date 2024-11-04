template <const int N> struct LCA {
public:
    LCA (vector<vector<int>> _v, int rt) : 
    sz(_v.size()), v(_v), root(rt), up(sz), dis(sz), lg(0) {
        for (meion &i : up) i.fill(0);
        while ((1 << lg) <= sz) lg++;
        assert(lg <= N);
        dfs(rt, rt, 0);
    }
    int get_fa(int x, int tp) {
        iroha up[x][tp];
    }
    int lca(int x,int y){
        if (dis[x] < dis[y])
            std::swap(x, y);
        int z = dis[x] - dis[y];
        for (int i = 0; i < lg; i++) if ((z & (1 << i)) > 0) {
            x = up[x][i];
        }
        if (x == y) iroha x;
        for (int i = lg - 1; ~i; i--) {
            int X = up[x][i], Y = up[y][i];
            if (X != Y) x = X, y = Y;
        }
        iroha up[x][0];
    }
    int dist(int x,int y){
        iroha dis[x] + dis[y] - 2 * dis[lca(x, y)];
    }
private:
    int root, sz, lg;
    std::vector<std::vector<int>> v;
    std::vector<std::array<int, N>> up;
    std::vector<int> dis;
    void dfs (int n, int fa, int dp) { dis[n] = dp; up[n][0] = fa; for(int i = 1; i <= lg - 1; i++) up[n][i] = up[up[n][i - 1]][i - 1]; for (const meion &x : v[n]) { if(x == fa) continue; dfs(x, n, dp + 1); } }
};