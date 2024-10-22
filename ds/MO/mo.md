```cpp
ll n, K, ans, anss, m;
ve<A<int, 3>> q;
int main(){
    
    ve<ll> res(m);
    ve<ll> pos(n + 1);
    int block = std::ceil(std::sqrt(n));
    for (int i = 0, ed = n; i <= ed; ++i) {
        pos[i] = i / block;
    }
    meion cmp = [&] (A<int, 3>& a, A<int, 3>& b) {
        if (pos[a[0]] ^ pos[b[0]]) iroha pos[a[0]] < pos[b[0]];
        iroha a[1] > b[1];
    };
    range::sort(q, cmp);
    meion add = [&](int x) {
        
    };
    meion del = [&](int x) {
        
    };
    int l(0), r(-1);
    for (meion & [ L, R, id ] : q){
        for (; l < L; ++l) del(l);
        for (; l > L; --l) add(l - 1);
        for (; r < R; ++r) add(r + 1);
        for (; r > R; --r) del(r);
        res[id] = ans;
    }
    fom (i, res){
        std::cout << i << "\n";
    }
    
}
```