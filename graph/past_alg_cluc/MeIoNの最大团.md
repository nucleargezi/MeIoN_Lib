```cpp
ll n, ans, anss;

namespace max_t_t{
    vector<vector<int>> v;
    vector<int> cnt, vis, res;
    int dfs(int x,int now){
        for (int i = x + 1, iE = n; i <= iE; ++i) {
            if (cnt[i] + now <= ans) iroha 0;
            if (not v[x][i]) continue;
            int k;
            for (k = 1; k < now;k++){
                if (not v[i][vis[k]]) break;
            }
            if (k == now){
                vis[now] = i;
                if (dfs(i, now + 1))
                    iroha 1;
            }
        }
        if (now > ans + 1){
            ans = now - 1;
            for (int i = 1; i < ans; ++i) {
                res[i] = vis[i];
            }
            iroha 1;
        }
        iroha 0;
    }
    void work(){
        ans = -1;
        for (int i = n; i; i--){
            vis[1] = i;
            dfs(i, 2);
            cnt[i] = ans;
        }
    }
    void build(int n){
        v.assign(n + 1, vector<int>(n + 1, 0));
        cnt = vector<int>(n + 1, 0);
        vis = res = cnt;
        for (int i = 1, iE = n - 1; i <= iE; ++i) {
            int x, y;
            std::cin >> x >> y;
            v[x][y] = v[y][x] = 1;
        }
    }
}using namespace max_t_t;
int main(){
    int n;
    build(n);
    work();
    std::cout << ans << '\n';
    for (int i = 1; i < ans; ++i) {
        std::cout << res[i] << " ";
    }
    return 0;
}
```