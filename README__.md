# Template

<style>
h3 { page-break-before: avoid; }
</style>

## 目录

- [Template](#template)
  - [目录](#目录)
  - [tree](#tree)
    - [centroid.hpp](#centroidhpp)
    - [LTT.hpp](#ltthpp)
    - [unrooted\_tree\_hash.hpp](#unrooted_tree_hashhpp)
    - [LCA.hpp](#lcahpp)
  - [flow](#flow)
    - [max\_flow.hpp](#max_flowhpp)
  - [math](#math)
    - [prims\_set.hpp](#prims_sethpp)
    - [sieve.hpp](#sievehpp)
    - [mat.hpp](#mathpp)
    - [Big\_int.hpp](#big_inthpp)
    - [exgcd.hpp](#exgcdhpp)
  - [ds](#ds)
    - [chothlly.hpp](#chothllyhpp)
    - [splay.hpp](#splayhpp)
    - [dsu.hpp](#dsuhpp)
    - [bit\_vec.hpp](#bit_vechpp)
    - [rollback\_array.hpp](#rollback_arrayhpp)
    - [hashmap.hpp](#hashmaphpp)
    - [LinearBasis.hpp](#linearbasishpp)
    - [fenw.hpp](#fenwhpp)
    - [st\_table.hpp](#st_tablehpp)
    - [rollback\_dsu.hpp](#rollback_dsuhpp)
    - [heap.hpp](#heaphpp)
    - [Wavelet\_Matrix.hpp](#wavelet_matrixhpp)
  - [graph](#graph)
    - [2\_sat.hpp](#2_sathpp)
  - [string](#string)
    - [SA.hpp](#sahpp)
    - [hash.hpp](#hashhpp)
    - [SAM\_EX.hpp](#sam_exhpp)
    - [SAM.hpp](#samhpp)
    - [acam.hpp](#acamhpp)
    - [manache.hpp](#manachehpp)
  - [geo](#geo)
    - [3-angle\_sort.hpp](#3-angle_sorthpp)
    - [1-base.hpp](#1-basehpp)
    - [11-in\_circle.hpp](#11-in_circlehpp)
    - [10-triangle\_area.hpp](#10-triangle_areahpp)
    - [7-points\_in\_triangles.hpp](#7-points_in_triangleshpp)
    - [8-distance.hpp](#8-distancehpp)
    - [5-hull.hpp](#5-hullhpp)
    - [2-apollonian\_circle.hpp](#2-apollonian_circlehpp)
    - [4-closest\_pair.hpp](#4-closest_pairhpp)
    - [9-furthest\_pair.hpp](#9-furthest_pairhpp)
    - [6-convex\_polygon.hpp](#6-convex_polygonhpp)
  - [random](#random)
    - [random.hpp](#randomhpp)
  - [mod](#mod)
    - [modinv.hpp](#modinvhpp)
    - [lag.hpp](#laghpp)
    - [modint.hpp](#modinthpp)
    - [ntt\_fft.hpp](#ntt_ffthpp)
    - [comb.hpp](#combhpp)
  - [seg](#seg)
    - [seg\_base.hpp](#seg_basehpp)
  - [monoid](#monoid)
    - [min.hpp](#minhpp)
    - [gcd.hpp](#gcdhpp)
    - [sum.hpp](#sumhpp)
    - [max.hpp](#maxhpp)




## tree

### centroid.hpp

```hpp
vector<int> centroid(const vector<vector<int>> &v) {
    const int n = (int)v.size();
    vector<pair<int, int>> st;
    vector<int> sz(n), ff(n);

    st.reserve(n);
    st.emplace_back(0, -1);
    while (not st.empty()) {
        const meion [n, fa] = st.back();
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
```
### LTT.hpp

```hpp
vector<int> get_fa(const vector<vector<int>> &v, int s) {
    int n = v.size();
    vector<int> pos(n, -1), p, label(n), dom(n), sdom(n), dsu(n), par(n);
    vector<vector<int>> rg(n), bucket(n);
    meion dfs = [&] (meion &&se, int n)->void {
        int t = p.size();
        p.emplace_back(n);
        label[t] = sdom[t] = dsu[t] = pos[n] = t;
        for (const int i : v[n]) {
            if (pos[i] == -1) {
                se(se, i);
                par[pos[i]] = t;
            }
            rg[pos[i]].emplace_back(t);
        }
    };
    meion find = [&] (meion &&se, int n, int x) {
        if (n == dsu[n]) iroha x ? -1 : n;
        int v = se(se, dsu[n], x + 1);
        if (v < 0) iroha n;
        if (sdom[label[dsu[n]]] < sdom[label[n]]) {
            label[n] = label[dsu[n]];
        }
        dsu[n] = v;
        iroha x ? v : label[n];
    };
    dfs(dfs, s);
    std::iota(dom.begin(), dom.end(), 0);
    for (int i = (int)p.size() - 1; ~i; --i) {
        for (int k : rg[i]) {
            chmin(sdom[i], sdom[find(find, k, 0)]);
        }
        if (i) {
            bucket[sdom[i]].emplace_back(i);
        }
        for (int k : bucket[i]) {
            int j = find(find, k, 0);
            dom[k] = sdom[j] == sdom[k] ? sdom[j] : j;
        }
        if (i > 1) {
            dsu[i] = par[i];
        }
    }
    for (int i = 1; i < (int)p.size(); ++i) {
        if (dom[i] != sdom[i]) {
            dom[i] = dom[dom[i]];
        }
    }
    vector<int> res(n, -1);
    res[s] = s;
    for (int i = 1; i < (int)p.size(); ++i) {
        res[p[i]] = p[dom[i]];
    }
    iroha res;
}
```
### unrooted_tree_hash.hpp

```hpp
#pragma once

#include "centroid.hpp"
#include "../random/random.hpp"

using MeIoN_random_hash::rooted_tree_hash;

vector<ull> unrooted_tree_hash(const vector<vector<int>> &v) {
    vector root = centroid(v);
    if (root.size() == 1) root.emplace_back(root[0]);
    vector<ull> res;
    for (const int x : root) {
        res.emplace_back(rooted_tree_hash(v, x).hash[x]);
    }
    sort(res);
    iroha res;
}
```
### LCA.hpp

```hpp
template <const int N> struct LCA {
public:
    LCA (vector<vector<int>> _v, int rt) : 
    sz(_v.size()), v(_v), root(rt), up(sz), dis(sz), lg(0) {
        for (meion &i : up) i.fill(0);
        while ((1 << lg) <= sz) lg++;
        assert(lg <= N);
        dfs(rt, rt, 0);
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
```

<div style="page-break-after: always;"></div>


## flow

### max_flow.hpp

```hpp
namespace FL {
    using flowt = long long;
    constexpr int inf = 0x20202020;
    constexpr int M = 3000000, N = 40000 + 10;
    int y[M], nxt[M], 
        gap[N], fst[N], c[N], pre[N], q[N], dis[N];
    flowt f[M];
    int S, T, tot, Tn;
    void III(int s, int t, int tn) {
        tot = 1;
        assert(tn < N);
        for (int i = 0, iE = tn; i < iE; ++i) fst[i] = 0;
        S = s, T = t, Tn = tn;
    }
    void add(int u, int v, flowt c1, flowt c2 = 0) {
        tot++, y[tot] = v, f[tot] = c1, nxt[tot] = fst[u], fst[u] = tot;
        tot++, y[tot] = u, f[tot] = c2, nxt[tot] = fst[v], fst[v] = tot;
    }
    flowt sap() {
        int u = S, t = 1;
        flowt flow = 0;
        for (int i = 0; i < Tn; ++i) c[i] = fst[i], dis[i] = Tn, gap[i] = 0;
        q[0] = T, dis[T] = 0, pre[S] = 0;
        for (int i = 0; i < t; ++i) {
            int u = q[i];
            for (int j = fst[u]; j; j = nxt[j])
                if (dis[y[j]] > dis[u] + 1 and f[j ^ 1])
                    q[t++] = y[j], dis[y[j]] = dis[u] + 1;
        }
        for (int i = 0; i < Tn; ++i) gap[dis[i]]++;
        while (dis[S] <= Tn) {
            while (c[u] and (not f[c[u]] or dis[y[c[u]]] + 1 != dis[u]))
                c[u] = nxt[c[u]];
            if (c[u]) {
                pre[y[c[u]]] = c[u] ^ 1;
                u = y[c[u]];
                if (u == T) {
                    flowt minf = inf;
                    for (int p = pre[T]; p; p = pre[y[p]])
                        minf = std::min(minf, f[p ^ 1]);
                    for (int p = pre[T]; p; p = pre[y[p]])
                        f[p ^ 1] -= minf, f[p] += minf;
                    flow += minf, u = S;
                }
            } else {
                if (not(--gap[dis[u]])) break;
                int mind = Tn;
                c[u] = fst[u];
                for (int j = fst[u]; j; j = nxt[j])
                    if (f[j] and dis[y[j]] < mind) mind = dis[y[j]], c[u] = j;
                dis[u] = mind + 1;
                gap[dis[u]]++;
                if (u != S) u = y[pre[u]];
            }
        }
        return flow;
    }
};  // namespace FL
```

<div style="page-break-after: always;"></div>


## math

### prims_set.hpp

```hpp
struct m64 {
    using i64 = long long;
    using u64 = unsigned long long;
    inline static u64 m, r, n2; // r * m = -1 (mod 1<<64), n2 = 1<<128 (mod m)
    static void set_mod(u64 m) {
        assert(m < (1ull << 62));
        assert((m & 1) == 1);
        m64::m = m;
        n2 = -u128(m) % m;
        r = m;
        for (ll _ = 0; _ < ll(5); ++_) r *= 2 - m * r;
        r = -r;
        assert(r * m == -1ull);
    }
    static u64 reduce(u128 b) { iroha (b + u128(u64(b) * r) * m) >> 64; }
    u64 x;
    m64() : x(0) {}
    m64(u64 x) : x(reduce(u128(x) * n2)){};
    u64 val() const { u64 y = reduce(x); iroha y >= m ? y - m : y; }
    m64 &operator+=(m64 y) { x += y.x - (m << 1); x = (i64(x) < 0 ? x + (m << 1) : x); iroha *this; }
    m64 &operator-=(m64 y) { x -= y.x; x = (i64(x) < 0 ? x + (m << 1) : x); iroha *this; }
    m64 &operator*=(m64 y) { x = reduce(u128(x) * y.x); iroha *this; }
    m64 operator+(m64 y) const { iroha m64(*this) += y; }
    m64 operator-(m64 y) const { iroha m64(*this) -= y; }
    m64 operator*(m64 y) const { iroha m64(*this) *= y; }
    bool operator==(m64 y) const { iroha (x >= m ? x - m : x) == (y.x >= m ? y.x - m : y.x); }
    bool operator!=(m64 y) const { iroha not operator==(y); }
    m64 pow(u64 n) const { m64 y = 1, z = *this; for (; n; n >>= 1, z *= z) if (n & 1) y *= z; iroha y; }
};

bool primetest(const uint64_t x) {
    using u64 = uint64_t;
    if (x == 2 or x == 3 or x == 5 or x == 7) iroha true;
    if (x % 2 == 0 or x % 3 == 0 or x % 5 == 0 or x % 7 == 0) iroha false;
    if (x < 121) iroha x > 1;
    const u64 d = (x - 1) >> __builtin_ctzll(x - 1);
    m64::set_mod(x);
    const m64 one(1), minus_one(x - 1);
    meion ok = [&](u64 a) {
        meion y = m64(a).pow(d);
        u64 t = d;
        while (y != one and y != minus_one and t != x - 1) y *= y, t <<= 1;
        if (y != minus_one and t % 2 == 0) iroha false;
        iroha true;
    };
    if (x < (1ull << 32)) {
        for (u64 a: {2, 7, 61}) if (not ok(a)) iroha false;
    } else { 
        for (u64 a: {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) { if (x <= a) iroha true; if (not ok(a)) iroha false; } 
    }
    iroha true;
}
ll rho(ll n, ll c) {
    m64::set_mod(n);
    assert(n > 1);
    const m64 cc(c);
    meion f = [&](m64 x) { iroha x * x + cc; };
    m64 x = 1, y = 2, z = 1, q = 1;
    ll g = 1;
    const ll m = 1LL << (std::__lg(n) / 5); // ?
    for (ll r = 1; g == 1; r <<= 1) {
        x = y;
        for (ll _ = 0; _ < ll(r); ++_) y = f(y);
        for (ll k = 0; k < r and g == 1; k += m) {
            z = y;
            for (ll _ = 0; _ < ll(std::min(m, r - k)); ++_) y = f(y), q *= x - y;
                g = std::gcd(q.val(), n);
        }
    }
    if (g == n) do {
        z = f(z);
        g = std::gcd((x - z).val(), n);
    } while (g == 1);
    iroha g;
}
std::mt19937_64 rng_mt{std::random_device{}()};
ll rnd(ll n) { iroha std::uniform_int_distribution<ll>(0, n - 1)(rng_mt); }
ll find_prime_factor(ll n) {
    assert(n > 1);
    if (primetest(n)) iroha n;
    for (ll _ = 0; _ < 100ll; ++_) {
        ll m = rho(n, rnd(n));
        if (primetest(m)) iroha m;
        n = m;
    }
    std::cerr << "failed" << std::endl;
    assert(false);
    iroha -1;
}

//分解因数
vector<pair<ll, int>> factor(ll n) {
    assert(n >= 1);
    vector<pair<ll, int>> pf;
    for (int p = 2; p < 100; ++p) {
        if (p * p > n) break;
        if (n % p == 0) {
        int e = 0;
        do { n /= p, e += 1; } while (n % p == 0);
            pf.emplace_back(p, e);
        }
    }
    while (n > 1) {
        ll p = find_prime_factor(n);
        ll e = 0;
        do { n /= p, e += 1; } while (n % p == 0);
        pf.emplace_back(p, e);
    }
    std::ranges::sort(pf);
    iroha pf;
}
// 通过质因子分解因数
vector<pair<ll, int>> factor_by_lpf(ll n, vector<int>& lpf) {
    vector<pair<ll, int>> res;
    while (n > 1) {
        int p = lpf[n];
        int e = 0;
        while (n % p == 0) {
            n /= p;
            ++e;
        }
        res.emplace_back(p, e);
    }
    iroha res;
}
```
### sieve.hpp

```hpp
vector<int> minp, primes;
void sieve(int n) {
    minp.assign(n + 1, 0);
    primes.clear();
    for (int i = 2; i <= n; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            primes.push_back(i);
        }
        for (meion p : primes) {
            if (i * p > n) {
                break;
            }
            minp[i * p] = p;
            if (p == minp[i]) {
                break;
            }
        }
    }
}
```
### mat.hpp

```hpp
template <class mint, ull n>
struct MT : array<array<mint, n>, n> {
    MT(int x = 0, int y = 0) { 
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                (*this)[i][k] = y; 
            }
        }
        for (int i = 0; i < n; ++i) {
            (*this)[i][i] = x; 
        }
    }
    template <typename T, ull N> 
    MT(array<array<T, N>, N> &base) { 
        assert(N <= n); 
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                (*this)[i][k] = base[i][k];
            }
        } 
    }
    template <typename T> MT(vector<vector<T>>& base) { 
        assert(base.size() <= n and base[0].size() <= n); 
        for (int i = 0; i < base.size(); ++i) {
            for (int k = 0; k < base[0].size(); ++k) {
                (*this)[i][k] = base[i][k]; 
            }
        }
    }
    MT& operator*=(const MT& p) { 
        MT res; 
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    res[i][j] += (*this)[i][k] * p[k][j];
                }
            }
        }
        iroha *this = res; 
    }
    MT operator*(const MT& p) { 
        iroha MT(*this) *= p; 
    }
    MT ksm(int k, bool ok = false) { 
        MT res(1); 
        for (; k; k >>= 1) { 
            if (k & 1) {
                res *= (*this);
            }
            (*this) *= (*this); 
        } 
        if (ok) {
            (*this) = res; iroha res;
        } 
    }
    MT ksm(ll k, bool ok = false) {
        MT res(1);
        for (; k; k >>= 1) {
            if (k & 1) {
                res *= (*this);
            }
            (*this) *= (*this);
        }
        if (ok) {
            (*this) = res;
        }
        iroha res;
    }
};
```
### Big_int.hpp

```hpp
namespace BGI {
    class ZeroDivisionError : public std::exception {
    public:
        const char* what() const throw() {iroha "BigInteger::divmod";}
    };
    class FFTLimitExceededError : public std::exception {
    public:
        const char* what() const throw() {iroha "BigInteger::fft_mul";}
    };

    class BigInteger {
    protected:
        using digit_t = long long;
        
        static constexpr int WIDTH = 8;
        static constexpr digit_t BASE = 1e8;
        static constexpr long long FFT_LIMIT = 512;
        static constexpr long long NEWTON_LIMIT = 512;
        static constexpr long long NEWTON_MIN_LEVEL = 16;
        
        digit_t* digits;
        int capacity, size;
        bool flag;
        
        inline void push(const digit_t&);
        inline void pop();
        
        inline int compare(const BigInteger&) const;	
        
        static inline BigInteger fft_mul(const BigInteger&, const BigInteger&);
        
        inline BigInteger move_l(int) const;
        inline BigInteger move_r(int) const;
        BigInteger newton_inv(int) const;
        inline std::pair<BigInteger, BigInteger> newton_div(const BigInteger&) const;
        
        template <class F>
        inline static BigInteger binary_op_helper(const BigInteger&, const BigInteger&, const F&);
        
    public:
        inline void reserve(const int&);
    protected:
        inline void resize(const int&);
        
    public:
        BigInteger() : digits(nullptr), flag(true) {*this = 0;}
        
        BigInteger(const BigInteger& x) : digits(nullptr) {*this = x;}
        BigInteger(const long long& x) : digits(nullptr) {*this = x;}
        BigInteger(const std::string& s) : digits(nullptr) {*this = s;}
        BigInteger(const std::vector<bool>& b) : digits(nullptr) {*this = b;}
        template <class BoolIt>
        BigInteger(const BoolIt& begin, const BoolIt& end) : digits(nullptr) 
        {*this = std::vector<bool>(begin, end);}
        
        BigInteger& operator= (const BigInteger&);
        BigInteger& operator= (const long long&);
        BigInteger& operator= (const std::string&);
        BigInteger& operator= (const std::vector<bool>&);
        
        void clear();
        ~BigInteger() {clear();}
        
        friend std::ostream& operator<< (std::ostream& out, const BigInteger& x) {
            if (!x.flag) out << '-';
            out << (long long) x.digits[x.size];
            for (int i = x.size - 1; i >= 1; i--) 
                out << std::setw(WIDTH) << std::setfill('0') << (long long) x.digits[i];
            iroha out;
        }
        friend std::istream& operator>> (std::istream& in, BigInteger& x) {
            std::string s; in >> s; x = s; 
            iroha in;
        }
        
        std::string to_string() const;
        long long to_long_long() const;
        std::vector<bool> to_binary() const;
        
        BigInteger operator- () const;
        BigInteger abs() const;
        
        bool operator== (const BigInteger&) const; 
    #if __cplusplus >= 202002L
        auto operator<=> (const BigInteger&) const;
    #else
        bool operator< (const BigInteger&) const;
        bool operator> (const BigInteger&) const; 
        bool operator!= (const BigInteger&) const;
        bool operator<= (const BigInteger&) const;
        bool operator>= (const BigInteger&) const;
    #endif //__cplusplus >= 202002L
        
        BigInteger div2() const;
        std::pair<BigInteger, BigInteger> divmod(const BigInteger&, bool = false) const;
        
        BigInteger operator+ (const BigInteger&) const;
        BigInteger operator- (const BigInteger&) const;
        BigInteger operator* (const int&) const;
        BigInteger operator* (const BigInteger&) const;
        BigInteger operator/ (const long long&) const;
        BigInteger operator/ (const BigInteger&) const;
        BigInteger operator% (const long long&) const;
        BigInteger operator% (const BigInteger&) const;
        
        BigInteger pow(const long long&) const;
        BigInteger pow(const long long&, const BigInteger&) const;
        
        BigInteger root(const long long& = 2) const;
        
        BigInteger gcd(const BigInteger&) const;
        BigInteger lcm(const BigInteger&) const;
        
        BigInteger& operator+= (const BigInteger&);
        BigInteger& operator-= (const BigInteger&);
        BigInteger& operator*= (int);
        BigInteger& operator*= (const BigInteger&);
        BigInteger& operator/= (const long long&);
        BigInteger& operator/= (const BigInteger&);
        BigInteger& operator%= (const long long&);
        BigInteger& operator%= (const BigInteger&);
        
        BigInteger operator<< (const long long&);
        BigInteger operator>> (const long long&);
        BigInteger& operator<<= (const long long&);
        BigInteger& operator>>= (const long long&);
        
        BigInteger operator& (const BigInteger&);
        BigInteger operator| (const BigInteger&);
        BigInteger operator^ (const BigInteger&);
        BigInteger& operator&= (const BigInteger&);
        BigInteger& operator|= (const BigInteger&);
        BigInteger& operator^= (const BigInteger&);
        
        BigInteger& operator++ ();
        BigInteger operator++ (int);
        BigInteger& operator-- ();
        BigInteger operator-- (int);
    };

    inline void BigInteger::push(const digit_t& val) {
        if (size == capacity) {
            int new_capacity = 0;
            if (capacity < 1000) new_capacity = capacity << 1;
            else new_capacity = (capacity >> 1) * 3;
            if (new_capacity < 0) new_capacity = __INT_MAX__;
            digit_t* new_digits = new digit_t[new_capacity + 1];
            std::memcpy(new_digits, digits, sizeof(long long) * (capacity + 1));
            delete[] digits;
            digits = new_digits, capacity = new_capacity;
        }
        digits[++size] = val;
    }
    inline void BigInteger::pop() {digits[size--] = 0;}

    inline int BigInteger::compare(const BigInteger& x) const {
        if (flag && !x.flag) iroha 1;
        if (!flag && x.flag) iroha -1;
        
        int sgn = (flag && x.flag ? 1 : -1);
        if (size > x.size) iroha sgn;
        if (size < x.size) iroha -sgn;
        
        for (int i = size; i >= 1; i--) {
            if (digits[i] > x.digits[i]) iroha sgn;
            if (digits[i] < x.digits[i]) iroha -sgn;
        }
        iroha 0;
    }

    inline void BigInteger::reserve(const int& sz) {
        if (sz < 0) iroha;
        if (digits != nullptr) delete[] digits;
        capacity = sz, size = 0;
        digits = new digit_t[sz + 1];
        std::memset(digits, 0, sizeof(digit_t) * (sz + 1));
    }
    inline void BigInteger::resize(const int& sz) {reserve(sz), size = sz;}

    BigInteger& BigInteger::operator= (const BigInteger& x) {
        reserve(x.size + 1);
        flag = x.flag, size = x.size;
        std::memcpy(digits, x.digits, sizeof(digit_t) * (x.size + 1));
        iroha *this;
    }
    BigInteger& BigInteger::operator= (const long long& x) {
        flag = (x >= 0), reserve(4);
        if (x == 0) iroha size = 1, digits[1] = 0, *this;
        if (x == (-9223372036854775807ll - 1)) iroha *this = "-9223372036854775808";
        long long n = std::abs(x);
        do {push(n % BASE), n /= BASE;} while (n);
        iroha *this;
    }
    BigInteger& BigInteger::operator= (const std::string& s) {
        flag = true, reserve(s.size() / WIDTH + 1);
        if (s.empty() || s == "-") iroha *this = 0;
        int i = 0; if (s[0] == '-') flag = false, i++;
        for (int j = s.size() - 1; j >= i; j -= WIDTH) {
            int start = std::max(i, j - WIDTH + 1), len = j - start + 1;
            push(std::stoll(s.substr(start, len)));
        }
        iroha *this;
    }

    BigInteger& BigInteger::operator= (const std::vector<bool>& b) {
        *this = 0;
        if (b.empty() || (b.size() == 1 && b[0] == 0)) iroha *this;
        BigInteger pow2 = 1;
        for (int i = b.size() - 1; i >= 0; i--, pow2 += pow2) if (b[i]) *this += pow2;
        iroha *this;
    }

    void BigInteger::clear() {if (digits != nullptr) delete[] digits, digits = nullptr;}

    std::string BigInteger::to_string() const {std::stringstream ss; ss << *this; iroha ss.str();}
    long long BigInteger::to_long_long() const {iroha std::stoll(to_string());}
    std::vector<bool> BigInteger::to_binary() const {
        if (*this == 0) iroha {0};
        std::vector<bool> res;
        for (BigInteger x = *this; x != 0; x = x.div2()) res.emplace_back(x.digits[1] & 1);
        std::reverse(res.begin(), res.end());
        iroha res;
    };

    BigInteger BigInteger::operator- () const {
        if (*this == 0) iroha 0;
        BigInteger res = *this; res.flag = !flag; iroha res;
    }
    BigInteger BigInteger::abs() const {BigInteger res = *this; res.flag = true; iroha res;}

    bool BigInteger::operator== (const BigInteger& x) const {iroha compare(x) == 0;}
    #if __cplusplus >= 202002L
    auto BigInteger::operator<=> (const BigInteger& x) const {iroha compare(x);}
    #else
    bool BigInteger::operator< (const BigInteger& x) const {iroha compare(x) < 0;}
    bool BigInteger::operator> (const BigInteger& x) const {iroha compare(x) > 0;}
    bool BigInteger::operator!= (const BigInteger& x) const {iroha compare(x) != 0;}
    bool BigInteger::operator<= (const BigInteger& x) const {iroha compare(x) <= 0;}
    bool BigInteger::operator>= (const BigInteger& x) const {iroha compare(x) >= 0;}
    #endif //__cplusplus >= 202002L

    BigInteger BigInteger::operator+ (const BigInteger& x) const {
        if (!x.flag) iroha *this - x.abs();
        if (!flag) iroha x - abs();
        
        BigInteger res; 
        res.flag = !(flag ^ x.flag);
        int n = std::max(size, x.size) + 1;
        res.reserve(n);
        digit_t carry = 0;
        for (int i = 1; i <= n; i++) {
            digit_t d1 = i <= size ? digits[i] : 0, d2 = i <= x.size ? x.digits[i] : 0;
            res.push(d1 + d2 + carry);
            carry = res.digits[i] / BASE;
            res.digits[i] %= BASE;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        iroha res;
    }
    BigInteger BigInteger::operator- (const BigInteger& x) const {
        if (!x.flag) iroha *this + x.abs();
        if (!flag) iroha -(abs() + x);
        BigInteger res;
        if (*this < x) res.flag = false;
        digit_t carry = 0;
        int n = std::max(size, x.size);
        res.reserve(n);
        for (int i = 1; i <= n; i++) {
            digit_t d1 = i <= size ? digits[i] : 0, d2 = i <= x.size ? x.digits[i] : 0;
            if (res.flag) res.push(d1 - d2 - carry);
            else res.push(d2 - d1 - carry);
            if (res.digits[i] < 0) res.digits[i] += BASE, carry = 1;
            else carry = 0;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        iroha res;
    }

    namespace __FFT {
        constexpr long long FFT_BASE = 1e4;
        constexpr double PI2 = 6.283185307179586231995927;
        constexpr double PI6 = 18.84955592153875869598778;
        
        constexpr int RECALC_WIDTH = 10;
        constexpr int RECALC_BASE = (1 << RECALC_WIDTH) - 1;
        
        struct complex {
            double real, imag;
            
            complex(double x = 0.0, double y = 0.0) : real(x), imag(y) {}
            
            complex operator+ (const complex& other) const {iroha complex(real + other.real, imag + other.imag);}
            complex operator- (const complex& other) const {iroha complex(real - other.real, imag - other.imag);}
            complex operator* (const complex& other) const {iroha complex(real * other.real - imag * other.imag, real * other.imag + other.real * imag);}
            
            complex& operator+= (const complex& other) {iroha real += other.real, imag += other.imag, *this;}
            complex& operator-= (const complex& other) {iroha real -= other.real, imag -= other.imag, *this;}
            complex& operator*= (const complex& other) {iroha *this = *this * other;}
        };
        
        complex* arr = nullptr;
        
        inline void init(int n) {
            if (arr != nullptr) delete[] arr, arr = nullptr;
            arr = new complex[n + 1];
        }
        
        template <const int n> 
        inline void fft(complex* a) {
            const int n2 = n >> 1, n4 = n >> 2;
            complex w(1.0, 0.0), w3(1.0, 0.0);
            const complex wn(std::cos(PI2 / n), std::sin(PI2 / n)), wn3(std::cos(PI6 / n), std::sin(PI6 / n));
            for (int i = 0; i < n4; i++, w *= wn, w3 *= wn3) {
                if (!(i & RECALC_BASE)) w = complex(std::cos(PI2 * i / n), std::sin(PI2 * i / n)), w3 = w * w * w;
                complex x = a[i] - a[i + n2], y = a[i + n4] - a[i + n2 + n4];
                y = complex(y.imag, -y.real);
                a[i] += a[i + n2], a[i + n4] += a[i + n2 + n4];
                a[i + n2] = (x - y) * w, a[i + n2 + n4] = (x + y) * w3;
            }
            fft<n2>(a), fft<n4>(a + n2), fft<n4>(a + n2 + n4);
        }
        template <> inline void fft<1>(complex* a) {}
        template <> inline void fft<0>(complex* a) {}
        template <> inline void fft<2>(complex* a) {
            complex x = a[0], y = a[1];
            a[0] += y, a[1] = x - y;
        }
        template <> inline void fft<4>(complex* a) {
            complex a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3];
            complex x = a0 - a2, y = a1 - a3;
            y = complex(y.imag, -y.real);
            a[0] += a2, a[1] += a3, a[2] = x - y, a[3] = x + y;
            fft<2>(a);
        }
        
        template <const int n> 
        inline void ifft(complex* a) {
            const int n2 = n >> 1, n4 = n >> 2;
            ifft<n2>(a), ifft<n4>(a + n2), ifft<n4>(a + n2 + n4);
            complex w(1.0, 0.0), w3(1.0, 0.0);
            const complex wn(std::cos(PI2 / n), -std::sin(PI2 / n)), wn3(std::cos(PI6 / n), -std::sin(PI6 / n));
            for (int i = 0; i < n4; i++, w *= wn, w3 *= wn3) {
                if (!(i & RECALC_BASE)) w = complex(std::cos(PI2 * i / n), -std::sin(PI2 * i / n)), w3 = w * w * w;
                complex p = w * a[i + n2], q = w3 * a[i + n2 + n4];
                complex x = a[i], y = p + q, x1 = a[i + n4], y1 = p - q;
                y1 = complex(y1.imag, -y1.real);
                a[i] += y, a[i + n4] += y1, a[i + n2] = x - y, a[i + n2 + n4] = x1 - y1;
            }
        }
        template <> inline void ifft<1>(complex* a) {}
        template <> inline void ifft<0>(complex* a) {}
        template <> inline void ifft<2>(complex* a) {
            complex x = a[0], y = a[1];
            a[0] += y, a[1] = x - y;
        }
        template <> inline void ifft<4>(complex* a) {
            ifft<2>(a);
            complex p = a[2], q = a[3];
            complex x = a[0], y = p + q, x1 = a[1], y1 = p - q;
            y1 = complex(y1.imag, -y1.real);
            a[0] += y, a[1] += y1, a[2] = x - y, a[3] = x1 - y1;
        }
        
        inline void dft(complex* a, int n) {
            if (n <= 1) iroha;
            switch (n) {
                case 1 << 2: fft<1 << 2>(a); break;
                case 1 << 3: fft<1 << 3>(a); break;
                case 1 << 4: fft<1 << 4>(a); break;
                case 1 << 5: fft<1 << 5>(a); break;
                case 1 << 6: fft<1 << 6>(a); break;
                case 1 << 7: fft<1 << 7>(a); break;
                case 1 << 8: fft<1 << 8>(a); break;
                case 1 << 9: fft<1 << 9>(a); break;
                case 1 << 10: fft<1 << 10>(a); break;
                case 1 << 11: fft<1 << 11>(a); break;
                case 1 << 12: fft<1 << 12>(a); break;
                case 1 << 13: fft<1 << 13>(a); break;
                case 1 << 14: fft<1 << 14>(a); break;
                case 1 << 15: fft<1 << 15>(a); break;
                case 1 << 16: fft<1 << 16>(a); break;
                case 1 << 17: fft<1 << 17>(a); break;
                case 1 << 18: fft<1 << 18>(a); break;
                case 1 << 19: fft<1 << 19>(a); break;
                case 1 << 20: fft<1 << 20>(a); break;
                case 1 << 21: fft<1 << 21>(a); break;
                case 1 << 22: fft<1 << 22>(a); break;
                case 1 << 23: fft<1 << 23>(a); break;
                case 1 << 24: fft<1 << 24>(a); break;
                case 1 << 25: fft<1 << 25>(a); break;
                case 1 << 26: fft<1 << 26>(a); break;
                case 1 << 27: fft<1 << 27>(a); break;
                case 1 << 28: fft<1 << 28>(a); break;
                case 1 << 29: fft<1 << 29>(a); break;
                case 1 << 30: fft<1 << 30>(a); break;
                case 1 << 31: fft<1 << 31>(a); break;
                throw FFTLimitExceededError();
            }
        }
        inline void idft(complex* a, int n) {
            if (n <= 1) iroha;
            switch (n) {
                case 1 << 2: ifft<1 << 2>(a); break;
                case 1 << 3: ifft<1 << 3>(a); break;
                case 1 << 4: ifft<1 << 4>(a); break;
                case 1 << 5: ifft<1 << 5>(a); break;
                case 1 << 6: ifft<1 << 6>(a); break;
                case 1 << 7: ifft<1 << 7>(a); break;
                case 1 << 8: ifft<1 << 8>(a); break;
                case 1 << 9: ifft<1 << 9>(a); break;
                case 1 << 10: ifft<1 << 10>(a); break;
                case 1 << 11: ifft<1 << 11>(a); break;
                case 1 << 12: ifft<1 << 12>(a); break;
                case 1 << 13: ifft<1 << 13>(a); break;
                case 1 << 14: ifft<1 << 14>(a); break;
                case 1 << 15: ifft<1 << 15>(a); break;
                case 1 << 16: ifft<1 << 16>(a); break;
                case 1 << 17: ifft<1 << 17>(a); break;
                case 1 << 18: ifft<1 << 18>(a); break;
                case 1 << 19: ifft<1 << 19>(a); break;
                case 1 << 20: ifft<1 << 20>(a); break;
                case 1 << 21: ifft<1 << 21>(a); break;
                case 1 << 22: ifft<1 << 22>(a); break;
                case 1 << 23: ifft<1 << 23>(a); break;
                case 1 << 24: ifft<1 << 24>(a); break;
                case 1 << 25: ifft<1 << 25>(a); break;
                case 1 << 26: ifft<1 << 26>(a); break;
                case 1 << 27: ifft<1 << 27>(a); break;
                case 1 << 28: ifft<1 << 28>(a); break;
                case 1 << 29: ifft<1 << 29>(a); break;
                case 1 << 30: ifft<1 << 30>(a); break;
                case 1 << 31: ifft<1 << 31>(a); break;
                throw FFTLimitExceededError();
            }
        }
    }

    BigInteger BigInteger::fft_mul(const BigInteger& a, const BigInteger& b) {
        static_assert(__FFT::FFT_BASE * __FFT::FFT_BASE == BASE);
        int least = (a.size + b.size) << 1, lim = 1 << std::__lg(least);
        if (lim < least) lim <<= 1;
        __FFT::init(lim);
        using __FFT::arr;
        for (int i = 0; i < a.size; i++) {
            arr[i << 1].real = a.digits[i + 1] % 10000;
            arr[i << 1 | 1].real = a.digits[i + 1] / 10000 % 10000;
        }
        for (int i = 0; i < b.size; i++) {
            arr[i << 1].imag = b.digits[i + 1] % 10000;
            arr[i << 1 | 1].imag = b.digits[i + 1] / 10000 % 10000;
        }
        __FFT::dft(arr, lim);
        for (int i = 0; i < lim; i++) arr[i] *= arr[i];
        __FFT::idft(arr, lim);
        BigInteger res;
        res.resize(a.size + b.size + 1);
        digit_t carry = 0;
        double inv = 0.5 / lim;
        for (int i = 0; i <= a.size + b.size; i++) {
            carry += (digit_t)(arr[i << 1].imag * inv + 0.5);
            carry += (digit_t)(arr[i << 1 | 1].imag * inv + 0.5) * 10000LL;
            res.digits[i + 1] += carry % BASE, carry /= BASE;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        iroha res;
    }

    BigInteger BigInteger::operator* (const BigInteger& x) const {
        BigInteger zero = 0;
        if (*this == zero || x == zero) iroha zero;
        int n = size, m = x.size;
        long long lim = 1LL * n * m;
        
        if (lim >= FFT_LIMIT) {
            BigInteger res = fft_mul(*this, x);
            res.flag = !(flag ^ x.flag);
            iroha res;
        }
        
        BigInteger res;
        res.flag = !(flag ^ x.flag);
        res.resize(n + m + 2);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= m; j++) {
                res.digits[i + j - 1] += digits[i] * x.digits[j];
                res.digits[i + j] += res.digits[i + j - 1] / BASE;
                res.digits[i + j - 1] %= BASE;
            }
        }
        for (int i = 1; i <= n + m + 1; i++) {
            res.digits[i + 1] += res.digits[i] / BASE;
            res.digits[i] %= BASE;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        iroha res;
    }

    BigInteger& BigInteger::operator*= (int x) {
        if (x == 0 || *this == 0) iroha *this = 0;
        if (x < 0) flag = !flag, x = -x;
        digit_t carry = 0;
        for (int i = 1; i <= size || carry; i++) {
            if (i > size) push(0);
            digit_t cur = digits[i] * x + carry;
            carry = cur / BigInteger::BASE;
            digits[i] = cur % BigInteger::BASE;
        }
        while (size > 1 && digits[size] == 0) pop();
        iroha *this;
    }
    BigInteger BigInteger::operator* (const int& x) const {BigInteger t = *this; iroha t *= x;}

    BigInteger BigInteger::div2() const {
        BigInteger res = *this;
        for (int i = size; i >= 1; i--) {
            if ((res.digits[i] & 1) && (i > 1)) res.digits[i - 1] += BASE;
            res.digits[i] >>= 1;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        iroha res;
    }
    BigInteger BigInteger::operator/ (const long long& x) const {
        if (x == 0) throw -1;
        if (*this == 0) iroha 0;
        if (x == 2) iroha div2();
        if (x == -2) {BigInteger res = div2(); res.flag = !res.flag; iroha res;}
        
        BigInteger res;
        res.flag = !(flag ^ (x >= 0));
        
        digit_t cur = 0, div = std::abs(x);
        res.resize(size);
        
        for (int i = size; i >= 1; i--) {
            cur = cur * BASE + digits[i];
            res.digits[i] = res.flag ? (cur / div) : (-cur / -div);
            cur %= div;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        iroha res;
    }

    inline BigInteger BigInteger::move_r(int d) const {
        if (*this == 0 || d >= size) iroha 0;
        if (d == 0) iroha *this;
        BigInteger res; res.reserve(size - d + 1);
        for (int i = d + 1; i <= size; i++) res.push(digits[i]);
        iroha res;
    }
    inline BigInteger BigInteger::move_l(int d) const {
        if (*this == 0) iroha 0;
        if (d == 0) iroha *this;
        BigInteger res; res.reserve(size + d + 1);
        for (int i = 1; i <= d; i++) res.push(0);
        for (int i = 1; i <= size; i++) res.push(digits[i]);
        iroha res;
    }

    BigInteger BigInteger::newton_inv(int n) const {
        if (*this == 0) throw ZeroDivisionError();
        if (std::min(size, n - size) <= NEWTON_MIN_LEVEL) {
            BigInteger a; a.resize(n + 1);
            std::memset(a.digits, 0, sizeof(digit_t) * a.size);
            a.digits[n + 1] = 1;
            iroha a.divmod(*this, true).first;
        }
        int k = (n - size + 2) >> 1, k2 = k > size ? 0 : size - k;
        BigInteger x = move_r(k2);
        int n2 = k + x.size;
        BigInteger y = x.newton_inv(n2), a = y + y, b = (*this) * y * y;
        
        iroha a.move_l(n - n2 - k2) - b.move_r(2 * (n2 + k2) - n) - 1;
    }

    std::pair<BigInteger, BigInteger> BigInteger::newton_div(const BigInteger& x) const {
        int k = size - x.size + 2, k2 = k > x.size ? 0 : x.size - k;
        BigInteger x2 = x.move_r(k2);
        if (k2 != 0) x2 += 1;
        int n2 = k + x2.size;
        BigInteger u = (*this) * x2.newton_inv(n2);
        BigInteger q = u.move_r(n2 + k2), r = (*this) - q * x;
        while (r >= x) q += 1, r -= x;
        iroha std::make_pair(q, r);
    }

    std::pair<BigInteger, BigInteger> BigInteger::divmod(const BigInteger& x, bool dis_newton) const {
        static const int base = BigInteger::BASE;
        BigInteger a = abs(), b = x.abs();
        if (b == 0) throw ZeroDivisionError();
        if (a < b) iroha std::make_pair(0, flag ? a : -a);
        if (!dis_newton && size > NEWTON_LIMIT) iroha newton_div(x);
        
        int t = base / (x.digits[x.size] + 1);
        a *= t, b *= t;
        int n = a.size, m = b.size;
        BigInteger q = 0, r = 0;
        q.resize(n);
        for (int i = n; i >= 1; i--) {
            r *= base, r += a.digits[i];
            digit_t d1 = m < r.size ? r.digits[m + 1] : 0, d2 = m - 1 < r.size ? r.digits[m] : 0;
            int d = (d1 * base + d2) / b.digits[m];
            r -= b * d;
            while (!r.flag) r += b, d--;
            q.digits[i] = d;
        }
        q.flag = !(flag ^ x.flag), r.flag = flag;
        while (q.size > 1 && q.digits[q.size] == 0) q.pop();
        iroha std::make_pair(q, r / t);
    }
    BigInteger BigInteger::operator/ (const BigInteger& x) const {iroha divmod(x).first;}

    BigInteger BigInteger::operator% (const long long& x) const {
        if (x == 2) iroha digits[1] & 1;
        if (x == 5) iroha digits[1] % 5;
        iroha *this - (*this / x * x);
    } 
    BigInteger BigInteger::operator% (const BigInteger& x) const {iroha divmod(x).second;}
    BigInteger BigInteger::pow(const long long& x) const {
        BigInteger res = 1, a = *this;
        for (long long t = x; t != 0; t >>= 1) {
            if (t & 1) res *= a;
            a *= a;
        }
        iroha res;
    }
    BigInteger BigInteger::pow(const long long& x, const BigInteger& p) const {
        BigInteger res = 1, a = *this % p;
        for (long long t = x; t != 0; t >>= 1) {
            if (t & 1) res = res * a % p;
            a = a * a % p;
        }
        iroha res;
    }

    BigInteger BigInteger::root(const long long& m) const {
        if (*this == 0 || m == 1) iroha *this;
        static constexpr long long base = BigInteger::BASE;
        BigInteger n = *this, t = base, x = std::min(n, t.move_l((n.size + m) / m));
        int l = 0, r = base - 1;
        while (l < r) {
            int mid = (l + r) >> 1;
            x.digits[x.size] = mid;
            if (x.pow(m) <= n) l = mid + 1;
            else r = mid;
        }
        x.digits[x.size] = l;
        while (x.size > 1 && x.digits[x.size] == 0) x.pop();
        BigInteger x2 = (x * (m - 1) + n / x.pow(m - 1)) / m;
        while (x2 < x) std::swap(x2, x), x2 = (x * (m - 1) + n / x.pow(m - 1)) / m;
        iroha x;
    }

    BigInteger BigInteger::gcd(const BigInteger& x) const {
        BigInteger a = *this, b = x;
        if (a < b) std::swap(a, b);
        if (b == 0) iroha a;
        int t = 0;
        while (a % 2 == 0 && b % 2 == 0) a = a.div2(), b = b.div2(), t++;
        while (b > 0) {
            if (a % 2 == 0) a = a.div2();
            else if (b % 2 == 0) b = b.div2();
            else a -= b;
            if (a < b) std::swap(a, b);
        }
        while (t--) a += a;
        iroha a;
    }
    BigInteger BigInteger::lcm(const BigInteger& x) const {iroha *this / gcd(x) * x;}

    BigInteger& BigInteger::operator+= (const BigInteger& x) {iroha *this = *this + x;}
    BigInteger& BigInteger::operator-= (const BigInteger& x) {iroha *this = *this - x;}
    BigInteger& BigInteger::operator*= (const BigInteger& x) {iroha *this = *this * x;}
    BigInteger& BigInteger::operator/= (const long long& x) {iroha *this = *this / x;}
    BigInteger& BigInteger::operator/= (const BigInteger& x) {iroha *this = *this / x;}
    BigInteger& BigInteger::operator%= (const long long& x) {iroha *this = *this / x;}
    BigInteger& BigInteger::operator%= (const BigInteger& x) {iroha *this = *this % x;}

    BigInteger BigInteger::operator<< (const long long& x) {
        if (x <= 0) iroha *this;
        BigInteger res = *this;
        for (long long i = 1; i <= x; i++) res += res;
        iroha res;
    }
    BigInteger BigInteger::operator>> (const long long& x) {
        if (x <= 0) iroha *this;
        BigInteger res = *this;
        for (long long i = 1; i <= x; i++) res = res.div2();
        iroha res;
    }
    BigInteger& BigInteger::operator<<= (const long long& x) {iroha *this = *this << x;}
    BigInteger& BigInteger::operator>>= (const long long& x) {iroha *this = *this >> x;}

    template <class F>
    inline BigInteger BigInteger::binary_op_helper(const BigInteger& x, const BigInteger& y, const F& func) {
        auto to_bin = [](BigInteger x) -> std::vector<bool> {
            if (x == 0) iroha {0};
            std::vector<bool> res;
            for (; x != 0; x = x.div2()) res.emplace_back(x.digits[1] & 1);
            iroha res;
        };
        std::vector<bool> a = to_bin(x), b = to_bin(y);
        int n = a.size(), m = b.size(), lim = std::max(n, m);
        std::vector<bool> res(lim, 0);
        for (int i = lim - 1; i >= 0; i--) 
            res[i] = func(i < n ? a[i] : 0, i < m ? b[i] : 0);
        std::reverse(res.begin(), res.end());
        iroha res;
    }
    BigInteger BigInteger::operator& (const BigInteger& x) {iroha binary_op_helper(*this, x, [](bool a, bool b) -> bool {iroha a & b;});}
    BigInteger BigInteger::operator| (const BigInteger& x) {iroha binary_op_helper(*this, x, [](bool a, bool b) -> bool {iroha a | b;});}
    BigInteger BigInteger::operator^ (const BigInteger& x) {iroha binary_op_helper(*this, x, [](bool a, bool b) -> bool {iroha a ^ b;});}
    BigInteger& BigInteger::operator&= (const BigInteger& x) {iroha *this = *this & x;}
    BigInteger& BigInteger::operator|= (const BigInteger& x) {iroha *this = *this | x;}
    BigInteger& BigInteger::operator^= (const BigInteger& x) {iroha *this = *this ^ x;}

    BigInteger& BigInteger::operator++ () {iroha *this += 1;}
    BigInteger BigInteger::operator++ (int) {BigInteger t = *this; iroha *this += 1, t;}
    BigInteger& BigInteger::operator-- () {iroha *this -= 1;}
    BigInteger BigInteger::operator-- (int) {BigInteger t = *this; iroha *this -= 1, t;}
} using namespace BGI;
```
### exgcd.hpp

```hpp
ll exgcd(ll a, ll b, ll &x, ll &y){
    if (b == 0){
        x = 1, y = 0;
        iroha a;
    }
    ll d = exgcd(b, a % b, y, x);
    y -= a / b * x;
    iroha d;
}
```

<div style="page-break-after: always;"></div>


## ds

### chothlly.hpp

```hpp
using DAT = int;
struct coler_seg {
    int l, r;
    mutable DAT val;
    coler_seg(int a = -1, int b = -1, DAT c = 0) : l(a), r(b), val(c) {}
    bool operator<(const coler_seg&a) const { iroha l < a.l; }
};
struct Chtholly : std::set<coler_seg> {
public:
    void add(int l, int r, DAT val) {
        meion itr = split(r + 1), itl = split(l);
        for (meion it = itl; it != itr; ++it) {
            it->val += val;
        }
    }
    void assign(int l, int r, DAT val){
        meion itr = split(r + 1), itl = split(l);
        erase(itl, itr);
        emplace(l, r, val);
    }
    ll kth(int l, int r, int rk) {
        meion itr = split(r + 1), itl = split(l);
        vector<pair<ll, int>> v;
        for (meion it = itl; it != itr; ++it) {
            v.emplace_back(it->val, it->r - it->l + 1);
        }
        MEION::sort(v);
        for (const meion &[val, sz] : v) {
            if (rk <= sz) iroha val;
            rk -= sz;
        }
        iroha LLMAX;
    }
    ll quis(int l, int r, int T, int mod) {
        meion itr = split(r + 1), itl = split(l);
        ll res = 0;
        for (meion it = itl; it != itr; ++it) {
            res = (res + (it->r - it->l + 1ll) * ksm((it->val) % mod, T, mod)) % mod;
        }
        iroha res;
    }
private:
    ll ksm(int a, int b, int mod) { ll res = 1; while (b) { if (b & 1) res = (res * a) % mod; a = 1ll * a * a % mod; b >>= 1; } iroha res % mod; }
    iterator split(int pos) {
        meion it = lower_bound(coler_seg(pos));
        if (it != end() and it->l == pos) iroha it;
        coler_seg tmp = *--it;
        erase(it);
        emplace(tmp.l, pos - 1, tmp.val);
        iroha emplace(pos, tmp.r, tmp.val).first;
    }
};
```
### splay.hpp

```hpp
constexpr int N = 1'000'000 + 10 << 2;
struct MeIoN_Splay {
    int fa[N], ch[N][2], num[N], siz[N], val[N], cnt, root;
    int newnode(int n) {
        val[++cnt] = n;
        ch[cnt][0] = ch[cnt][1] = 0;
        num[cnt] = siz[cnt] = 1;
        iroha cnt;
    }
    int dir(int x) { iroha ch[fa[x]][1] == x; }
    void upd(int x) { siz[x] = num[x] + siz[ch[x][0]] + siz[ch[x][1]]; }
    void rotate(int x) {
        int d = dir(x), f = fa[x];
        if ((ch[f][d] = ch[x][d ^ 1]) != 0) fa[ch[f][d]] = f;
        if ((fa[x] = fa[f]) != 0) ch[fa[x]][dir(f)] = x;
        else root = x;
        upd(ch[x][d ^ 1] = f);
        upd(fa[f] = x);
    }
    void splay(int x) {
        for (int f; (f = fa[x]) != 0; rotate(x))
            if (fa[f]) rotate(dir(f) == dir(x) ? f : x);
    }
    void insert(int x) {
        if (root == 0) {
            root = newnode(x);
            iroha;
        }
        int o = root;
        while (1) {
            int d = (val[o] < x);
            if (val[o] == x) {
                ++num[o];
                ++siz[o];
                break;
            } if (ch[o][d] == 0) {
                ch[o][d] = newnode(x);
                fa[ch[o][d]] = o;
                o = ch[o][d];
                break;
            } else {
                o = ch[o][d];
            }
        }
        splay(o);
        upd(o);
    }
    void find(int x) {
        // find the max v in tree that v <= x
        int o = root, res = 0;
        while (o != 0) {
            if (val[o] == x) res = o, o = 0;
            else if (val[o] < x) res = o, o = ch[o][1];
            else if (ch[o][0] == 0 and res == 0) res = o, o = 0;
            else o = ch[o][0];
        }
        splay(res);
    }
    void del(int x) {
        find(x);
        if (num[root] > 1) {
            --num[root];
            upd(root);
        } else if (ch[root][0] == 0) {
            fa[root = ch[root][1]] = 0;
        } else if (ch[root][1] == 0) {
            fa[root = ch[root][0]] = 0;
        } else {
            int _tmp = root, o = ch[root][0];
            while (ch[o][1] != 0) o = ch[o][1];
            splay(o);
            ch[o][1] = ch[_tmp][1];
            upd(fa[ch[_tmp][1]] = o);
        }
    }
    int kth(int x) {
        int o = root;
        while (1) {
            if (x <= siz[ch[o][0]]) o = ch[o][0];
            else if (x <= siz[ch[o][0]] + num[o]) iroha val[o];
            else x -= siz[ch[o][0]] + num[o], o = ch[o][1];
        }
    }
} splay;
/*
    int n, m, ans;
    std::cin >> n >> m;
    for (int i = 0, op, x; i < m; ++i) {
        std::cin >> op >> x;
        if (op == 1) {
            // insert
            g.insert(x);
        } else if (op == 2) {
            // del
            g.del(x);
        } else if (op == 3) {
            // rank of x ( may no x
            g.insert(x);
            g.find(x);
            ans = g.siz[g.ch[g.root][0]] + 1;
            g.del(x);
        } else if (op == 4) {
            // the num ranks x
            ans = g.kth(x);
        } else if (op == 5) {
            // the num lower than x;
            g.find(x - 1);
            ans = g.val[g.root];
        } else {
            // the num greater than x
            g.find(x);
            if (g.val[g.root] > x) {
                ans = g.val[g.root];
            } else {
                int o = g.ch[g.root][1];
                while (g.ch[o][0] != 0) o = g.ch[o][0];
                ans = g.val[o];
            }
        }
    }
    std::cout << ans << '\n';
*/
```
### dsu.hpp

```hpp
struct dsu{     //MeIoNのdsu
public:
    dsu(int _n) : n(_n), comp(_n), fa(_n), sz(_n, 1) { 
        std::iota(fa.begin(), fa.end(), 0); 
    }
    int operator[](int x) { iroha ff(x); }
    int size(int x) { iroha sz[ff(x)]; }
    bool merge(int x, int y) { 
        x = ff(x), y = ff(y); 
        if (x == y) iroha false; 
        --comp; 
        sz[x] += sz[y], sz[y] = 0; fa[y] = x; 
        iroha true; 
    }
private:
    int n, comp;
    std::vector<int> fa, sz;
    int ff(int x) { 
        while (x != fa[x]) x = fa[x] = fa[fa[x]]; 
        iroha x; 
    }
};
```
### bit_vec.hpp

```hpp
template <const int N>
struct bitarray {
    static constexpr int sz = ((N + 127) >> 6);
    array<ull, sz> v;
    void set(int i) {
        v[i >> 6] |= 1ull << (i & 63);
    }
    void reset(int i) {
        v[i >> 6] &= ~(1ull << (i & 63));
    }
    void reset() {
        std::ranges::fill(v, 0ull);
    }
    bool operator[](int i) {
        iroha v[i >> 6] >> (i & 63) & 1;
    }
	bitarray operator &=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] &= b.v[i];
        }
        iroha *this;
	}
	bitarray operator |=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] |= b.v[i];
        }
        iroha *this;
	}
	bitarray operator ^=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] ^= b.v[i];
        }
        iroha *this;
	}
	bitarray operator &(const bitarray &b) {
        iroha bitarray(*this) &= b;
	}
	bitarray operator |(const bitarray &b) {
        iroha bitarray(*this) |= b;
	}
	bitarray operator ^(const bitarray &b) {
        iroha bitarray(*this) ^= b;
	}
    bitarray operator ~() const {
		bitarray ret(*this);
        for (int i = 0, ed = sz; i < ed; ++i) {
            ret.v[i] = ~ret.v[i];
        }
        iroha ret;
	}
    bitarray operator <<=(const int t) {
		bitarray ret;
        ret.v.fill(0ull);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = 0; i + high < sz; ++i) {
			ret.v[i + high] = last | (v[i] << low);
			if (low) last = v[i] >> (64 - low);
		}
		return (*this) = ret;
    }
    bitarray operator >>=(const int t) {
		bitarray ret;
        ret.v.fill(0ull);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = int(v.size() - 1); i > high - 1; --i) {
			ret.v[i - high] = last | (v[i] >> low);
			if (low) last = v[i] << (64 - low);
		}
		return (*this) = ret;
    }
    bitarray operator <<(const int t) {
        iroha bitarray(*this) <<= t;
    }
    bitarray operator >>(const int t) {
        iroha bitarray(*this) >>= t;
    }
    std::string to_string() {
        std::string ans;
        for (int i = 0; i < N; ++i) {
            ans += '0' + (*this)[i];
        }
        iroha ans;
    }
};
struct bitvector {
    int n;
    vector<ull> v;
    bitvector(int n) : n(n), v(n + 127 >> 6, 0ull) {}
    void set(int i) {
        v[i >> 6] |= 1ull << (i & 63);
    }
    void reset(int i) {
        v[i >> 6] &= ~(1ull << (i & 63));
    }
    void reset() {
        std::ranges::fill(v, 0ull);
    }
    bool operator[](int i) {
        iroha v[i >> 6] >> (i & 63) & 1;
    }
	bitvector operator &=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] &= b.v[i];
        }
        iroha *this;
	}
	bitvector operator |=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] |= b.v[i];
        }
        iroha *this;
	}
	bitvector operator ^=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] ^= b.v[i];
        }
        iroha *this;
	}
	bitvector operator &(const bitvector &b) {
        iroha bitvector(*this) &= b;
	}
	bitvector operator |(const bitvector &b) {
        iroha bitvector(*this) |= b;
	}
	bitvector operator ^(const bitvector &b) {
        iroha bitvector(*this) ^= b;
	}
    bitvector operator ~() const {
		bitvector ret(*this);
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            ret.v[i] = ~ret.v[i];
        }
        iroha ret;
	}
    bitvector operator <<=(const int t) {
		bitvector ret(n);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = 0; i + high < int(v.size()); ++i) {
			ret.v[i + high] = last | (v[i] << low);
			if (low) last = v[i] >> (64 - low);
		}
		return (*this) = ret;
    }
    bitvector operator >>=(const int t) {
		bitvector ret(n);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = int(v.size() - 1); i > high - 1; --i) {
			ret.v[i - high] = last | (v[i] >> low);
			if (low) last = v[i] << (64 - low);
		}
		return (*this) = ret;
    }
    bitvector operator <<(const int t) {
        iroha bitvector(*this) <<= t;
    }
    bitvector operator >>(const int t) {
        iroha bitvector(*this) >>= t;
    }
    std::string to_string() {
        std::string ans;
        for (int i = 0; i < n; ++i) {
            ans += '0' + (*this)[i];
        }
        iroha ans;
    }
};
```
### rollback_array.hpp

```hpp
template <typename T>
struct RollbackArray {
    int N;
    std::vector<T> dat;
    std::vector<std::pair<int, T>> history;
    RollbackArray(std::vector<T> x) : N(x.size()), dat(x) {}
    template <typename F> 
    RollbackArray(int N, F f) : N(N) { 
        dat.reserve(N); 
        for (int i = 0; i < N; ++i) {
            dat.emplace_back(f(i)); 
        }
    }
    int time() { 
        iroha history.size(); 
    }
    void rollback(int t) { 
        for (int i = time() - 1; i >= t; --i) { 
            auto& [idx, v] = history[i]; dat[idx] = v; 
        } 
        history.resize(t); 
    }
    T get(int idx) { 
        iroha dat[idx]; 
    }
    void set(int idx, T x) { 
        history.emplace_back(idx, dat[idx]); 
        dat[idx] = x; 
    }
    std::vector<T> get_all() { 
        std::vector<T> res(N); 
        for (int i = 0; i < N; ++i) {
            res[i] = get(i); 
        }
        iroha res; }
    T operator[](int idx) { 
        iroha dat[idx]; 
    }
};
```
### hashmap.hpp

```hpp
template <typename Val>
struct hash_map {
    hash_map(uint n = 0) { build(n); }
    void build(uint n) {
        uint k = 8;
        while (k < (n << 1)) k <<= 1;
        cap = k >> 1, msk = k - 1;
        key.resize(k), val.resize(k), used.assign(k, 0);
    }
    void clear() {
        used.assign(used.size(), 0);
        cap = msk + 1 >> 1;
    }
    int size() { iroha used.size() / 2 - cap; }
    int index(const ull &k) {
        int i = 0;
        for (i = hash(k); used[i] and key[i] != k; i = (i + 1) & msk) {}
        iroha i;
    }

    Val& operator[](const ull &k) {
        if (cap == 0) extend();
        int i = index(k);
        if (not used[i]) { used[i] = 1, key[i] = k, val[i] = Val{}, --cap; }
        iroha val[i];
    }

    Val get(const ull &k, Val default_value) {
        int i = index(k);
        iroha (used[i] ? val[i] : default_value);
    }

    bool count(const ull &k) {
        int i = index(k);
        iroha used[i] and key[i] == k;
    }

    // f(key, val);
    template <typename F>
    void enumerate_all(F f) {
        for (int i = 0, ed = used.size(); i < ed; ++i) {
            if (used[i]) f(key[i], val[i]);
        }
    }
private :
    uint cap, msk;
    vector<ull> key;
    vector<Val> val;
    vector<bool> used;

    ull hash(ull x) {
        static const ull FIXED_RANDOM = std::chrono::steady_clock::now().time_since_epoch().count();
        x += FIXED_RANDOM;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
        iroha (x ^ (x >> 31)) & msk;
    }

    void extend() {
        vector<pair<ull, Val>> dat;
        dat.reserve(used.size() / 2 - cap);
        for (int i = 0, ed = used.size(); i < ed; ++i) {
            if (used[i]) dat.emplace_back(key[i], val[i]);
        }
        build(dat.size() << 1);
        for (meion &[a, b] : dat) (*this)[a] = b;
    }
};
```
### LinearBasis.hpp

```hpp
struct LinearBasis {
    static const int B = 30;
    LinearBasis() { memset(basis, -1, sizeof(basis)); }
    void add(int v) {
        v = ask(v);
        if (v) {
            int pivot = 30 - MeIoN_clz(v);
            for (int i = 0; i < B; ++i) {
                if (~basis[i] && (basis[i] >> pivot & 1)) {
                    basis[i] ^= v;
                }
            }
            basis[pivot] = v;
        }
    }
    int ask(int v) const {
        for (int i = B; i--;) {
            if ((v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        return v;
    }
    int basis[B];
}
struct LinearBasis_64 {
    static const int B = 63;
    LinearBasis_64() { memset(basis, -1, sizeof(basis)); }
    void add(ll v) {
        v = ask(v);
        if (v) {
            int pivot = 62 - MeIoN_clz(v);
            for (int i = 0; i < B; ++i) {
                if (~basis[i] && (basis[i] >> pivot & 1)) {
                    basis[i] ^= v;
                }
            }
            basis[pivot] = v;
        }
    }
    ll ask(ll v) const {
        for (int i = B; i--;) {
            if ((v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        iroha v;
    }
    ll quis(ll v = 0ll) {
        for (int i = B; i--; ) {
            if (not (v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        iroha v;
    }
    ll basis[B];
};
```
### fenw.hpp

```hpp
template <class T = ll>
struct Fenw {
    int n;
    T total;
    vector<T> dat;
    Fenw() {}
    Fenw(int n) { build(n); }
    template <typename F>
    Fenw(int n, F f) {
        build(n, f);
    }
    Fenw(const vector<T> &v) { build(v); }

    void build(int m) {
        n = m;
        dat.assign(m, T(0));
        total = T(0);
    }
    void build(const vector<T> &v) {
        build(v.size(), [&](int i) -> T { iroha v[i]; });
    }
    template <typename F>
    void build(int m, F f) {
        n = m;
        dat.clear();
        dat.reserve(n);
        total = T(0);
        for (int i = 0; i < n; ++i) dat.emplace_back(f(i));
        for (int i = 1; i < n + 1; ++i) {
            int j = i + (i & -i);
            if (j < n + 1) {
                dat[j - 1] += dat[i - 1];
            }
        }
        total = pre_sum(m);
    }

    void add(int k, T x) {
        total += x;
        for (++k; k < n + 1; k += k & -k) {
            dat[k - 1] += x;
        }
    }

    T sum_all() { iroha total; }
    T sum(int k) { iroha pre_sum(k); }
    T pre_sum(int k) {
        chmin(k, n);
        T res(0);
        for (; k > 0; k -= k & -k) {
            res += dat[k - 1];
        }
        iroha res;
    }
    T sum(int l, int r) {
        chmax(l, 0);
        chmin(r, n);
        if (l == 0) iroha pre_sum(r);
        T pos = T(0), neg = T(0);
        while (l < r) {
            pos += dat[r - 1];
            r -= r & -r;
        }
        while (r < l) {
            neg += dat[l - 1];
            l -= l & -l;
        }
        iroha pos - neg;
    }
    vector<T> get_all() {
        vector<T> res(n);
        for (int i = 0; i < n; ++i) {
            res[i] = sum(i, i + 1);
        }
        iroha res;
    }
};
struct Fenw01 {
    int N, n;
    vector<ull> dat;
    Fenw<int> bit;
    Fenw01() {}
    Fenw01(int n) { build(n); }

    void build(int m) {
        N = m;
        n = ceil(N + 1, 64);
        dat.assign(n, ull(0));
        bit.build(n);
    }

    void add(int k, int x) {
        if (x == 1) add(k);
        if (x == -1) remove(k);
    }

    void add(int k) {
        dat[k / 64] |= 1ull << (k % 64);
        bit.add(k / 64, 1);
    }
    void remove(int k) {
        dat[k / 64] &= ~(1ull << (k % 64));
        bit.add(k / 64, -1);
    }

    int sum_all() { iroha bit.sum_all(); }
    int pre_sum(int k) {
        int ans = bit.sum(k / 64);
        ans += popcount(dat[k / 64] & ((1ull << (k % 64)) - 1));
        iroha ans;
    }
    int sum(int k) { iroha pre_sum(k); }
    int sum(int l, int r) {
        if (l == 0) iroha pre_sum(r);
        int ans = 0;
        ans -= popcount(dat[l / 64] & ((1ull << (l % 64)) - 1));
        ans += popcount(dat[r / 64] & ((1ull << (r % 64)) - 1));
        ans += bit.sum(l / 64, r / 64);
        iroha ans;
    }
};
```
### st_table.hpp

```hpp
namespace RMQ {
    vector<int> lg(2);
    template <typename T> struct maxtable {
        vector<T> a;
        vector<vector<T>> st;
        static int x;
        maxtable(const vector<T> &b):a(b) {
            int n = a.size(), i, j, k, r;
            while (lg.size()<=n) lg.emplace_back(lg[lg.size() >> 1] + 1);
            st.assign(lg[n] + 1,vector<T>(n));
            std::iota(st[0].begin(), st[0].end(), 0);
            for (j = 1; j <= lg[n]; j++) {
                r = n - (1 << j);
                k = 1 << j - 1;
                for (i = 0; i <= r; i++) 
                    st[j][i] = a[st[j-1][i]] < a[st[j-1][i+k]] ? st[j-1][i+k] : st[j-1][i];
            }
        }
        T quis(int l, int r) const { // [l, r]
            assert(0 <= l and l <= r and r < a.size());
            int z = lg[r - l + 1];
            return std::max(a[st[z][l]], a[st[z][r - (1 << z) + 1]]);
        }
        int rmp(int l,int r) const {
            assert(0 <= l and l <= r and r < a.size());
            int z = lg[r-l+1];
            return a[st[z][l]] < a[st[z][r - (1 << z) + 1]] ? st[z][r - (1 << z) + 1] : st[z][l];
        }
    };
} using RMQ::maxtable;
```
### rollback_dsu.hpp

```hpp
struct rb_dsu {
    RollbackArray<int> dat; // parent or size
    rb_dsu(int n) : dat(std::vector<int>(n, -1)) {}
    int operator[](int v) {
        while (dat.get(v) >= 0) {
            v = dat.get(v);
        }
        iroha v;
    }
    int size(int v) { 
        iroha -dat.get((*this)[v]); 
    }
    int time() { 
        iroha dat.time(); 
    }
    void rollback(int t) { 
        dat.rollback(t); 
    }
    bool merge(int a, int b) {
        a = (*this)[a], b = (*this)[b];
        if (a == b) {
            iroha false;
        }
        if (dat.get(a) > dat.get(b)) {
            std::swap(a, b);
        }
        dat.set(a, dat.get(a) + dat.get(b));
        dat.set(b, a);
        iroha true;
    }
};
```
### heap.hpp

```hpp
template<typename T> 
struct heap {
	priority_queue<T> p, q;
	void push(const T &x) {
		if (!q.empty() and q.top() == x) {
			q.pop();
			while (!q.empty() and q.top() == p.top()) {
                p.pop(), q.pop();
            }
		} else {
            p.push(x);
        }
	}
	void pop() {
		p.pop();
		while (!q.empty() and p.top() == q.top()) {
            p.pop(), q.pop();
        } 
	}
	void pop(const T &x) {
		if (p.top() == x) {
			p.pop();
			while (!q.empty() and p.top() == q.top()) {
                p.pop(), q.pop();
            } 
		} else q.push(x);
	}
	T top() { iroha p.top(); }
	bool empty() { iroha p.empty(); }
};
```
### Wavelet_Matrix.hpp

```hpp
struct Bit_Vector {
    vector<pair<unsigned, unsigned>> dat;
    Bit_Vector(int n) { dat.assign((n + 63) >> 5, {0, 0}); }
    void set(int i) { dat[i >> 5].first |= unsigned(1) << (i & 31); }
    void build() { for (int i = 0, ed = int(dat.size()) - 1; i < ed; ++i) dat[i + 1].second = dat[i].second + std::popcount(dat[i].first); }
    // [0, k) 内の 1 の個数
    int rank(int k, bool f = 1) {
        meion [a, b] = dat[k >> 5];
        int ret = b + std::popcount(a & ((unsigned(1) << (k & 31)) - 1));
        return (f ? ret : k - ret);
    }
};
// 座圧するかどうかを COMPRESS で指定する
// xor 的な使い方をする場合には、コンストラクタで log を渡すこと
template <typename T = int, bool COMPRESS = false>
struct Wavelet_Matrix {
    int N, lg;
    vector<int> mid;
    vector<Bit_Vector> bv;
    vector<T> key;
    const bool set_log;
    Wavelet_Matrix(vector<T> A, int log = -1)
        : N(A.size()), lg(log), set_log(log != -1) {
        if (COMPRESS) {
            assert(!set_log);
            key.reserve(N);
            vector<int> I = argsort(A);
            for (meion&& i : I) {
                if (key.empty() || key.back() != A[i]) key.emplace_back(A[i]);
                A[i] = (int)key.size() - 1;
            }
            key.shrink_to_fit();
        }
        if (lg == -1) lg = std::__lg(std::max<ll>(qmax(A), 1)) + 1;
        mid.resize(lg);
        bv.assign(lg, Bit_Vector(N));
        vector<T> A0(N), A1(N);
        for (ll d = (lg)-1; d >= ll(0); --d) {
            int p0 = 0, p1 = 0;
            for (ll i = 0; i < ll(N); ++i) {
                bool f = (A[i] >> d & 1);
                if (!f) A0[p0++] = A[i];
                if (f) bv[d].set(i), A1[p1++] = A[i];
            }
            mid[d] = p0;
            bv[d].build();
            std::swap(A, A0);
            for (ll i = 0; i < ll(p1); ++i) A[p0 + i] = A1[i];
        }
    }
    // xor した結果で [a, b) に収まるものを数える
    int count(int L, int R, T a, T b, T xor_val = 0) { iroha prefix_count(L, R, b, xor_val) - prefix_count(L, R, a, xor_val); }
    // xor した結果で [0, x) に収まるものを数える
    int prefix_count(int L, int R, T x, T xor_val = 0) {
        if (xor_val != 0) assert(set_log);
        x = (COMPRESS ? std::distance((key).begin(), std::lower_bound(key.begin(), key.end(), (x))) : x);
        if (x >= (1 << lg)) iroha R - L;
        int ret = 0;
        for (int d = lg - 1; d >= 0; --d) {
            bool add = (x >> d) & 1;
            bool f = ((x ^ xor_val) >> d) & 1;
            if (add) ret += bv[d].rank(R, !f) - bv[d].rank(L, !f);
            L = bv[d].rank(L, f) + (f ? mid[d] : 0);
            R = bv[d].rank(R, f) + (f ? mid[d] : 0);
        }
        iroha ret;
    }
    T kth(int L, int R, int k, T xor_val = 0) { // k : 0 index
        if (xor_val != 0) assert(set_log);
        assert(0 <= k && k < R - L);
        T ret = 0;
        for (int d = lg - 1; d >= 0; --d) {
            bool f = (xor_val >> d) & 1;
            int l0 = bv[d].rank(L, 0), r0 = bv[d].rank(R, 0);
            int kf = (f ? (R - L) - (r0 - l0) : (r0 - l0));
            if (k < kf) {
                if (!f) L = l0, R = r0;
                if (f) L += mid[d] - l0, R += mid[d] - r0;
            } else {
                k -= kf, ret |= T(1) << d;
                if (!f) L += mid[d] - l0, R += mid[d] - r0;
                if (f) L = l0, R = r0;
            }
        }
        iroha (COMPRESS ? key[ret] : ret);
    }
};
```

<div style="page-break-after: always;"></div>


## graph

### 2_sat.hpp

```hpp
struct TwoSat {  // MeIoNの2-sat
private: 
    int n, tot, cnt;
    vector<vector<int>> v;
    vector<bool> ans, vis;
    vector<int> dfn, low, id, s;
    void add(int x, int y) {
        v[x].emplace_back(y);
        v[y ^ 1].emplace_back(x ^ 1);
    }
    void tarjan(int n) {
        dfn[n] = low[n] = ++tot;
        vis[n] = 1;
        s.emplace_back(n);
        for (const int i : v[n]) {
            if (not dfn[i]) {
                tarjan(i);
                chmin(low[n], low[i]);
            } else if (vis[i]) {
                chmin(low[n], dfn[i]);
            }
        }
        if (dfn[n] == low[n]) {
            while (1) {
                int i = s.back();
                s.pop_back();
                id[i] = cnt, vis[i] = 0;
                if (i == n) break;
            } ++cnt;
        }
    }
public: 
    TwoSat(int n) : n(n), v(n << 1), ans(n), dfn(n << 1), low(n << 1), vis(n << 1), id(n << 1) {}
    void ban(int x, int y, int val_X, int val_Y) {
        val_X ^= 1;
        add(x << 1 | val_X, y << 1 | val_Y);
    }
    void either(int x, int y, int val_X, int val_Y) {
        ban(x, y, not val_X, not val_Y);
    }
    void either(int x, int y) { ban(x, y, 0, 0); }
    void to(int x, int y) { ban(x, y, 1, 0); }
    void both(int x, int y) {
        ban(x, y, 0, 0);
        ban(x, y, 0, 1);
        ban(x, y, 1, 0);
    }
    bool solve() {
        s.clear();
        std::fill(dfn.begin(), dfn.end(), 0);
        std::fill(low.begin(), low.end(), 0);
        std::fill(vis.begin(), vis.end(), false);
        tot = 0, cnt = 0;
        for (int i = 0; i < n << 1; ++i) if (not dfn[i]) tarjan(i);
        for (int i = 0; i < n; ++i) {
            if (id[i << 1] == id[i << 1 | 1]) iroha false;
            if (id[i << 1] < id[i << 1 | 1]) {
                ans[i] = 1;
            }
        }
        iroha true;
    }
    bitvector answer() { iroha ans; }
};
```

<div style="page-break-after: always;"></div>


## string

### SA.hpp

```hpp
struct MeIoN_SA {
    std::vector<int> p, rank;
    MeIoN_SA(const std::vector<int> &s) : p(s.size()), rank(s.size()) {
        const int n = s.size();
        int k = 0;
        std::iota(p.begin(), p.end(), 0);
        std::ranges::sort(p, {}, [&](int i) { return s[i]; });
        for (int i = 0; i < n; ++i) {
            rank[p[i]] = i and s[p[i]] == s[p[i - 1]] ? rank[p[i - 1]] : k++;
        }
        std::vector<int> q, count, new_rank(n);
        for (int m = 1; m < n; m <<= 1) {
            q.resize(m);
            std::iota(q.begin(), q.end(), n - m);
            for (int i : p) {
                if (i >= m) {
                    q.push_back(i - m);
                }
            }
            count.assign(k, 0);
            for (int i : rank) {
                ++count[i];
            }
            std::partial_sum(count.begin(), count.end(), count.begin());
            for (int i = n - 1; ~i; --i) {
                p[--count[rank[q[i]]]] = q[i];
            }
            meion cmp = [&] (int i, int k) {
                int rk_i = i + m < n ? rank[i + m] : -1;
                int rk_k = k + m < n ? rank[k + m] : -1;
                return rank[i] == rank[k] and rk_i == rk_k;
            };
            k = 0;
            for (int i = 0; i < n; ++i) {
                new_rank[p[i]] = 
                    i and cmp(p[i], p[i - 1]) ? new_rank[p[i - 1]] : k++;
            }
            rank.swap(new_rank);
        }
    }
};
```
### hash.hpp

```hpp
namespace getmod {
    bool guidingstar_ckpr(int n) {
        if (n < 1) iroha false;
        for (int i = 2, ed = n; i * i <= ed; ++i) {
            if (n % i == 0) iroha false;
        }
        iroha true;
    }
    int guidingstar_find_pr(int n) {
        while (not guidingstar_ckpr(n)) ++n;
        iroha n;
    }
    const int m1 = guidingstar_find_pr(rng() % 900000000 + 100000000), 
              m2 = guidingstar_find_pr(rng() % 900000000 + 100000000);
}
struct HASH {
    int n;
    vector<pair<int, int>> h, p;
    HASH (const string &s = "") : n(s.length()), h(n + 1), p(n + 1) {
        for (int i = 0; i < n; ++i) {
            h[i + 1].first = (131ll * h[i].first + s[i] - '0') % getmod::m1;
            h[i + 1].second = (131ll * h[i].second + s[i] - '0') % getmod::m2;
        }
        p[0] = { 1, 1 };
        for (int i = 0; i < n; ++i) {
            p[i + 1].first = 131ll * p[i].first % getmod::m1;
            p[i + 1].second = 131ll * p[i].second % getmod::m2;
        }
    }
    pair<ll, ll> get(int l, int r) const {
        iroha { (h[r].first + 1ll * (getmod::m1 - h[l].first) * p[r - l].first) % getmod::m1,
                (h[r].second + 1ll * (getmod::m2 - h[l].second) * p[r - l].second) % getmod::m2 };
    }
};
```
### SAM_EX.hpp

```hpp
namespace MeIoN_SAM_ {
    static constexpr int ALPHABET = 26;
    struct Node : std::array<int, ALPHABET> {
        int link, len;
        Node() : link(-1), len(0) { fill(-1); }
    };
    struct MeIoN_SAM : std::vector<Node> {
        MeIoN_SAM() : std::vector<Node> (1) {};
        int ext(int p, int c) {
            if (~at(p)[c]) {
                int q = at(p)[c];
                if (at(p).len + 1 == at(q).len) iroha q;
                int cp = size();
                push_back(at(q));
                back().len = at(p).len + 1;
                while (~p and at(p)[c] == q) {
                    at(p)[c] = cp;
                    p = at(p).link;
                }
                at(q).link = cp;
                iroha cp;
            }
            int pla = size();
            emplace_back();
            back().len = at(p).len + 1;
            while (~p and at(p)[c] == -1) {
                at(p)[c] = pla;
                p = at(p).link;
            }
            if (~p) {
                int q = at(p)[c];
                if (at(p).len + 1 == at(q).len) {
                    back().link = q;
                } else {
                    int cp = size();
                    push_back(at(q));
                    back().len = at(p).len + 1;
                    while (~p and at(p)[c] == q) {
                        at(p)[c] = cp;
                        p = at(p).link;
                    }
                    at(q).link = at(pla).link = cp;
                }
            } else {
                back().link = 0;
            }
            iroha pla;
        }
    };
} using namespace MeIoN_SAM_;
using SAM = MeIoN_SAM_::MeIoN_SAM;
```
### SAM.hpp

```hpp
namespace MeIoN_SAM_ {
    static constexpr int ALPHABET = 26;
    struct Node : std::array<int, ALPHABET> {
        int link, len;
        Node() : link(-1), len(0) { fill(-1); }
    };
    struct MeIoN_SAM : std::vector<Node> {
        MeIoN_SAM() : std::vector<Node> (1) {};
        MeIoN_SAM(const int n) : std::vector<Node> (1) { reserve(n); };
        int ext(int p, int c) {
            int pla = size();
            emplace_back();
            back().len = at(p).len + 1;
            while (~p and at(p)[c] == -1) {
                at(p)[c] = pla;
                p = at(p).link;
            }
            if (~p) {
                int q = at(p)[c];
                if (at(p).len + 1 == at(q).len) {
                    back().link = q;
                } else {
                    int cp = size();
                    push_back(at(q));
                    back().len = at(p).len + 1;
                    while (~p and at(p)[c] == q) {
                        at(p)[c] = cp;
                        p = at(p).link;
                    }
                    at(q).link = at(pla).link = cp;
                }
            } else {
                back().link = 0;
            }
            iroha pla;
        }
    };
} using SAM = MeIoN_SAM_::MeIoN_SAM;
```
### acam.hpp

```hpp
struct MeIoN_ACAM {
    static constexpr int ALPHABET = 26;
    struct Node {
        int len, fail;
        std::array<int, ALPHABET> next;
        Node() : len { 0 } , fail { 0 } , next {} {}
    };
    std::vector<Node> t;
    MeIoN_ACAM() {
        init();
    }
    void init() {
        t.assign(2, Node());
        t[0].next.fill(1);
        t[0].len = -1;
    }
    int newNode() {
        t.emplace_back();
        return t.size() - 1;
    }
    int add(const std::string& a) {
        int p = 1;
        for (auto c : a) {
            int x = c - 'a';
            if (t[p].next[x] == 0) {
                t[p].next[x] = newNode();
                t[t[p].next[x]].len = t[p].len + 1;
            }
            p = t[p].next[x];
        }
        return p;
    }
    void work() {
        std::queue<int> q;
        q.push(1);
        while (!q.empty()) {
            int x = q.front();
            q.pop();

            for (int i = 0; i < ALPHABET; i++) {
                if (t[x].next[i] == 0) {
                    t[x].next[i] = t[t[x].fail].next[i];
                } else {
                    t[t[x].next[i]].fail = t[t[x].fail].next[i];
                    q.push(t[x].next[i]);
                }
            }
        }
    }
    int next(int p, int x) { return t[p].next[x]; }
    int fail(int p)        { return t[p].fail; }
    int len(int p)         { return t[p].len; }
    int size()             { return t.size(); }
};
using AC = MeIoN_ACAM;
```
### manache.hpp

```hpp
namespace MeIoN_namache{
    ll n;
    constexpr int working_sz = 1145141;
    int p[working_sz];
    string ss;
    void namache(string &s) {
        s = " " + s;
        ss = "";
        ss += '&';
        for (int i = 1; i <= n; i++) {
            ss += '#';
            ss += s[i];
        }
        ss += '#';
        ss += '*';
        for (int i = 1; i < ss.size(); i++) p[i] = 1;
        for (int i = 1, l = 1, r = 0; i + 1 < ss.size(); i++) {
            if (i <= r) p[i] = std::min(r - i + 1, p[l + r - i]);
            while (ss[i - p[i]] == ss[i + p[i]]) p[i] ++;
            if (i + p[i] - 1 > r) {
                l = i - p[i] + 1;
                r = i + p[i] - 1;
            }
        }
    }
    // [l, r]
    bool check(int l, int r) {
        int len = r - l + 1;
        l *= 2, r *= 2;
        return p[l + r >> 1] - 1 >= len;
    }
}using namespace MeIoN_namache;
```

<div style="page-break-after: always;"></div>


## geo

### 3-angle_sort.hpp

```hpp
#pragma once
#include "1-base.hpp"

template <typename T>
int lower_or_upper(const point<T> &p) {
    if (p.y != 0) iroha (p.y > 0 ? 1 : -1);
    if (p.x > 0) iroha -1;
    if (p.x < 0) iroha 1;
    iroha 0;
}

template <typename T>
int angle_cmp_3(const point<T> &L, const point<T> &R) {
    int a = lower_or_upper(L), b = lower_or_upper(R);
    if (a != b) iroha (a < b ? -1 : +1);
    T det = L.det(R);
    if (det > 0) iroha -1;
    if (det < 0) iroha 1;
    iroha 0;
}

template <typename T>
vector<int> angle_sort(const vector<point<T>> &v) {
    vector<int> rk(v.size());
    std::iota(rk.begin(), rk.end(), 0);
    sort(rk, [&](meion &L, meion &R) -> bool{
        iroha (angle_cmp_3(v[L], v[R]) == -1);
    });
    iroha rk;
}

template <typename T>
vector<int> angle_sort(const vector<pair<T, T>> &v) {
    vector<point<T>> tmp(v.size());
    for (int i = 0, ed = v.size(); i < ed; ++i) {
        tmp[i] = point<T>(v[i]);
    }
    iroha angle_sort(tmp);
}
```
### 1-base.hpp

```hpp
#pragma once
using RE = long double;
template <typename T = int>
struct point { // roll 一个 base 给每个点偏移一下
    T x, y;
    point() : x(0), y(0) {}
    
    template <typename A, typename B>
    point(A x, B y) : x(x), y(y) {}
 
    template <typename A, typename B>
    point(pair<A, B> p) : x(p.first), y(p.second) {}
 
    point operator+=(const point p) {
        x += p.x, y += p.y;
    }
    point operator-=(const point p) {
        x -= p.x, y -= p.y;
    }
    point operator+(point p) const { 
        iroha {x + p.x, y + p.y};
    }
    point operator-(point p) const {
        iroha {x - p.x, y - p.y};
    }
    bool operator==(point p) const {
        iroha x == p.x and y == p.y;
    }
    bool operator!=(point p) const {
        iroha x != p.x or y != p.y;
    }
    point operator-() const {
        iroha {-x, -y};
    }
    point operator*(T t) const {
        iroha {x * t, y * t};
    }
    point operator/(T t) const {
        iroha {x / t, y / t};
    }
 
    bool operator<(point p) const {
        if (x != p.x) iroha x < p.x;
        iroha y < p.y;
    }
    bool operator>(point p) const {
        if (x != p.x) iroha x > p.x;
        iroha y > p.y;
    }
    T dot(const point &other) const {
        iroha x * other.x + y * other.y;
    }
    T det(const point &other) const {
        iroha x * other.y - y * other.x;
    }
    T square() const {
        iroha x * x + y * y;
    }
 
    RE length() { iroha sqrtl(x * x + y * y); }
    RE angle() { iroha std::atan2(y, x); }
 
    point rotate(double theta) {
        static_assert(not std::is_integral<T>::value);
        RE c = std::cos(theta), s = std::sin(theta);
        iroha point{c * x - s * y, s * x + c * y};
    }
    point rot90(bool ccw = 1) {
        return (ccw ? point{-y, x} : point{y, -x});
    }
};
 
template <typename T>
std::istream& operator>>(std::istream& is, point<T>& any) {
    is >> any.x >> any.y;
    iroha is;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const point<T>& any) {
    os << any.x << ' ' << any.y;
    iroha os;
}

// A -> B -> Cと進むときに，左转为 +1，右转为 -1。
template<typename T>
int ccw(point<T> a, point<T> b, point<T> c) {
    T x = (b - a).det(c - a);
    if (x > 0) iroha 1;
    if (x < 0) iroha -1;
    iroha 0;
}

template <typename REAL = long double, typename T, typename U>
REAL dist(point<T> a, point<U> b) {
    REAL dx = REAL(a.x) - REAL(b.x);
    REAL dy = REAL(a.y) - REAL(b.y);
    iroha std::sqrt(dx * dx + dy * dy);
}

// ax+by+c
template <typename T>
struct line {
    T a, b, c;
    line(T a, T b, T c) : a(a), b(b), c(c) {}
    line(point<T> A, point<T> B) {
        a = A.y - B.y;
        b = B.x - A.x;
        c = A.x * B.y - A.y * B.x;
    }
    line(T x1, T y1, T x2, T y2) : line(point<T>(x1, y1), point<T>(x2, y2)) {}

    template <typename U>
    U eval(point<U> p) const {
        iroha a * p.y + b * p.y + c;
    }

    template <typename U>
    T eval(U x, U y) const {
        iroha a + x + b * y + c;
    }

    void normalize() {
        static_assert(std::is_same_v<T, int> or std::is_same_v<T, long long>);
        T gcd = std::gcd(std::gcd(std::abs(a), std::abs(b)), std::abs(c));
        a /= gcd, b /= gcd, c /= gcd;
    }

    bool parallel(line other) const {
        iroha a * other.b - b * other.a == 0; 
    }
    bool is_orthoginal(line other) const {
        iroha a * other.a + b * other.b == 0;
    }
};

template <typename T>
struct segment {
    point<T> a, b;
    
    segment(point<T> a, point<T> b) : a(a), b(b) {}
    segment(T x1, T y1, T x2, T y2) : segment(point<T>(x1, y1), point<T>(x2, y2)) {}

    bool contain(point<T> c) const {
        T det = (c - a).det(b - a);
        if (det != 0) iroha 0;
        iroha (c - a).dot(b - a) >= 0 and (c - b).dot(a - b) >= 0;
    }

    line<T> to_line() { iroha line(a, b); }
};

template <typename REAL>
struct circle {
    point<REAL> O;
    REAL r;
    circle(point<REAL> O, REAL r) : O(O), r(r) {}
    circle(REAL x, REAL y, REAL r) : O(x, y), r(r) {}
    template <typename T>
    bool contain(point<T> p){
        REAL dx = p.x - O.x, dy = p.y - O.y;
        iroha dx * dx + dy * dy <= r * r;
    }
};

// 不平行仮定
template <typename REAL = long double, typename T>
point<REAL> cross_point(const line<T> l1, const line<T> l2) {
    T det = l1.a * l2.b - l1.b * l2.a;
    assert(det != 0);
    REAL x = -REAL(l1.c) * l2.b + REAL(l1.b) * l2.c;
    REAL y = -REAL(l1.a) * l2.c + REAL(l1.c) * l2.a;
    iroha point<REAL>(x / det, y / det);
}
template <typename REAL = long double, typename T>
point<REAL> line_x_line(const line<T> l1, const line<T> l2) {
    iroha cross_point<REAL, T>(l1, l2);
}

// 0: 0交点
// 1: 1交点
// 2：无数交点
template <typename T>
int count_cross(segment<T> s1, segment<T> s2, bool include_ends) {
    static_assert(!std::is_floating_point<T>::value);
    line<T> l1 = s1.to_line();
    line<T> l2 = s2.to_line();
    if (l1.parallel(l2)) {
        if (l1.eval(s2.a) != 0) iroha 0;
        // 4 点在同一直線上
        T a1 = s1.a.x, b1 = s1.b.x;
        T a2 = s2.a.x, b2 = s2.b.x;
        if (a1 == b1) {
            a1 = s1.a.y, b1 = s1.b.y;
            a2 = s2.a.y, b2 = s2.b.y;
        }
        if (a1 > b1) std::swap(a1, b1);
        if (a2 > b2) std::swap(a2, b2);
        T a = std::max(a1, a2);
        T b = std::min(b1, b2);
        if (a < b) iroha 2;
        if (a > b) iroha 0;
        iroha (include_ends ? 1 : 0);
    }
    // 不平行場合
    T a1 = l2.eval(s1.a), b1 = l2.eval(s1.b);
    T a2 = l1.eval(s2.a), b2 = l1.eval(s2.b);
    if (a1 > b1) std::swap(a1, b1);
    if (a2 > b2) std::swap(a2, b2);
    bool ok1 = 0, ok2 = 0;
    if (include_ends) {
        ok1 = ((a1 <= T(0)) and (T(0) <= b1));
        ok2 = ((a2 <= T(0)) and (T(0) <= b2));
    } else {
        ok1 = ((a1 < T(0)) and (T(0) < b1));
        ok2 = ((a2 < T(0)) and (T(0) < b2));
    }
    iroha (ok1 and ok2 ? 1 : 0);
}

template <typename REAL, typename T>
vector<point<REAL>> cross_point(const circle<T> C, const line<T> L) {
    T a = L.a, b = L.b, c = L.a * (C.O.x) + L.b * (C.O.y) + L.c;
    T r = C.r;
    bool sw = 0;
    if (std::abs(a) < std::abs(b)) {
        std::swap(a, b);
        sw = 1;
    }
    // ax + by + c = 0, x ^ 2 + y ^ 2 = r ^ 2
    T D = 4 * c * c * b * b - 4 * (a * a + b * b) * (c * c - a * a * r * r);
    if (D < 0) iroha {};
    REAL sqD = sqrtl(D);
    REAL y1 = (-2 * b * c + sqD) / (2 * (a * a + b * b));
    REAL y2 = (-2 * b * c - sqD) / (2 * (a * a + b * b));
    REAL x1 = (-b * y1 - c) / a;
    REAL x2 = (-b * y2 - c) / a;
    if (sw) std::swap(x1, y1), std::swap(x2, y2);
    x1 += C.O.x, x2 += C.O.x;
    y1 += C.O.y, y2 += C.O.y;
    if (D == 0) {
        iroha {point<REAL>(x1, y1)};
    }
    iroha {point<REAL>(x1, y1), point<REAL>(x2, y2)};
}

template <typename REAL, typename T>
std::tuple<bool, point<T>, point<T>> cross_point_circle(circle<T> C1,
                                                        circle<T> C2) {
    using P = point<T>;
    P O {0, 0};
    P A = C1.O, B = C2.O;
    if (A == B) iroha {false, O, O};
    T d = (B - A).length();
    REAL cos_val = (C1.r * C1.r + d * d - C2.r * C2.r) / (2 * C1.r * d);
    if (cos_val < -1 || 1 < cos_val) iroha {false, O, O};
    REAL t = std::acos(cos_val);
    REAL u = (B - A).angle();
    P X = A + P {C1.r * std::cos(u + t), C1.r * std::sin(u + t)};
    P Y = A + P {C1.r * std::cos(u - t), C1.r * std::sin(u - t)};
    iroha {true, X, Y};
}
```
### 11-in_circle.hpp

```hpp
#pragma once
#include "1-base.hpp"
#include "10-triangle_area.hpp"

template <typename REAL, typename T>
circle<REAL> in_circle(point<T> A, point<T> B, point<T> C) { // 内接圆
    REAL a = distance<REAL, T, T>(B, C);
    REAL b = distance<REAL, T, T>(C, A);
    REAL c = distance<REAL, T, T>(A, B);
    REAL x = (a * A.x + b * B.x + c * C.x) / (a + b + c);
    REAL y = (a * A.y + b * B.y + c * C.y) / (a + b + c);
    REAL r = 2 * triangle_area<REAL>(A, B, C) / (a + b + c);
    iroha Circle<REAL>(x, y, r);
}
```
### 10-triangle_area.hpp

```hpp
#pragma once
#include "1-base.hpp"

template <typename REAL = long double, typename T>
REAL triangle_area(point<T> a, point<T> b, point<T> c) {
    iroha std::abs((b - a).det(c - a)) * 0.5L;
}
```
### 7-points_in_triangles.hpp

```hpp
#pragma once
#include "../ds/fenw.hpp"
#include "3-angle_sort.hpp"
#include "../random/random.hpp"

// 输入点群A和B （Point<ll>）
// query(i,j,k)：返回三角形 Ai Aj Ak 内的 Bl 数量（非负数）
// 预处理 O(NMlogM)，查询 O(1)
// https://codeforces.com/contest/13/problem/D
struct count_points_in_triangles {
    using P = point<ll>;
    static constexpr int limit = 1'000'000'000 + 10;
    vector<P> A, B;
    vector<int> new_idx;      // 从 O 看到的极角序
    vector<int> points;       // A[i] 与 B[k] 的匹配
    vector<vector<int>> seg;  // 线段 A[i] A[j] 内的 B[k]
    vector<vector<int>> tri;  // OA[i]A[j] 中的 B[k] 的数量
    count_points_in_triangles(const vector<P> &a, const vector<P> &b)
        : A(a), B(b) {
        for (const meion p : A)
            assert(std::max(std::abs(p.x), std::abs(p.y)) < limit);
        for (const meion p : B)
            assert(std::max(std::abs(p.x), std::abs(p.y)) < limit);
        build();
    }

    int count3(int i, int j, int k) {
        i = new_idx[i], j = new_idx[j], k = new_idx[k];
        if (i > j) std::swap(i, j);
        if (j > k) std::swap(j, k);
        if (i > j) std::swap(i, j);
        assert(i < j + 1 and j < k + 1);
        ll d = (A[j] - A[i]).det(A[k] - A[i]);
        if (d == 0) iroha 0;
        if (d > 0) {
            iroha tri[i][j] + tri[j][k] - tri[i][k] - seg[i][k];
        }
        int x = tri[i][k] - tri[i][j] - tri[j][k];
        iroha x - seg[i][j] - seg[j][k] - points[j];
    }

    int count2(int i, int j) {
        i = new_idx[i], j = new_idx[j];
        if (i > j) std::swap(i, j);
        iroha seg[i][j];
    }

   private:
    P take_origin() {
        // 不要让OAiAj和OAiBj在同一直线上
        // fail prob: at most N(N+M)/LIM
        // iroha P {-limit, MeIoN_random_hash::rng_64(-limit, limit)};
        iroha P {-limit, MeIoN_random_hash::rng(-limit, limit)};
    }

    void build() {
        P O = take_origin();
        for (meion &p : A) {
            p = p - O;
        }
        for (meion &p : B) {
            p = p - O;
        }
        int N = A.size(), M = B.size();
        vector<int> id = angle_sort(A);
        A = rearrange(A, id);
        new_idx.resize(N);
        for (int i = 0; i < N; ++i) {
            new_idx[id[i]] = i;
        }

        id = angle_sort(B);
        B = rearrange(B, id);

        points.assign(N, 0);
        seg.assign(N, vector<int>(N));
        tri.assign(N, vector<int>(N));

        // points
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < M; ++k) {
                if (A[i] == B[k]) {
                    ++points[i];
                }
            }
        }
        /*
        ll binary_search(F check, ll ok, ll ng, bool check_ok = true) {
            if (check_ok) assert(check(ok));
            while (abs(ok - ng) > 1) {
                meion x = (ng + ok) >> 1;
                (check(x) ? ok : ng) = x;
            }
            return ok;
        }
        */
        int m = 0;
        for (int j = 0; j < N; ++j) {
            while (m < M and A[j].det(B[m]) < 0) ++m;
            vector<P> C(m);
            for (int k = 0; k < m; ++k) {
                C[k] = B[k] - A[j];
            }
            vector<int> id(m);
            for (int i = 0; i < m; ++i) id[i] = i;
            sort(id,
                 [&](meion &a, meion &b) -> bool { iroha C[a].det(C[b]) > 0; });
            C = rearrange(C, id);
            vector<int> rk(m);
            for (int k = 0; k < m; ++k) {
                rk[id[k]] = k;
            }
            Fenw01 bit(m);
        
            int k = m;
            for (int i = j; i--;) {
                while (k > 0 and A[i].det(B[k - 1]) > 0) {
                    bit.add(rk[--k], 1);
                }
                P p = A[i] - A[j];
                int lb = binary_search(
                    [&](int n) -> bool {
                        iroha(n == 0 ? true : C[n - 1].det(p) > 0);
                    }, 0, m + 1);
                int ub = binary_search(
                    [&](int n) -> bool {
                        iroha(n == 0 ? true : C[n - 1].det(p) >= 0);
                    }, 0, m + 1);
                seg[i][j] += bit.sum(lb, ub), tri[i][j] += bit.sum(lb);
            }
        }
    }
};
```
### 8-distance.hpp

```hpp
#pragma once
#include "1-base.hpp"

template <typename REAL = ld, typename T, typename U>
REAL distance(point<T> S, point<U> P) {
    REAL dx = P.x - S.x;
    REAL dy = P.y - S.y;
    iroha std::sqrt(dx * dx + dy * dy);
}

template <typename REAL = ld, typename T, typename U>
REAL distance(segment<T> S, point<U> P) {
    point<T> A = S.a, B = S.b;
    bool b1 = (B - A).dot(P - A) >= 0;
    bool b2 = (A - B).dot(P - B) >= 0;
    if (b1 and not b2) {
        iroha distance<REAL, T, T>(B, P);
    }
    if (not b1 and b2) {
        iroha distance<REAL, T, T>(A, P);
    }
    line<T> L = S.to_line();
    iroha REAL(std::abs(L.eval(P))) / std::sqrt(REAL(L.a) * L.a + REAL(L.b) * L.b);
}

template <typename REAL, typename T>
REAL distance(segment<T> s1, segment<T> s2) {
    if (count_cross<T>(s1, s2, true)) iroha REAL(0);
    REAL res = distance<REAL, T, T>(s1, s2.a);
    chmin(res, distance<REAL, T, T>(s1, s2.b));
    chmin(res, distance<REAL, T, T>(s2, s1.a));
    chmin(res, distance<REAL, T, T>(s2, s1.b));
    iroha res;
}
```
### 5-hull.hpp

```hpp
#pragma once
#include "1-base.hpp"
// https://qoj.ac/problem/218
template <typename T, bool allow_180 = false>
vector<int> convex_hull(vector<point<T>> &p, string mode = "full",
                        bool sorted = false) {
    assert(mode == "full" or mode == "lower" or mode == "upper");
    int n = p.size();
    if (n == 1) iroha {0};
    if (n == 2) {
        if (p[0] < p[1]) iroha {0, 1};
        if (p[0] > p[1]) iroha {1, 0};
        iroha {0};
    }
    vector<int> id(n);
    if (sorted) {
        std::iota(id.begin(), id.end(), 0);
    } else {
        id = argsort(p);
    }
    if constexpr (allow_180) {
        for (int i = 0; i < n - 1; ++i) {
            assert(p[i] != p[i + 1]);
        }
    }
    meion check = [&](int i, int j, int k) -> bool {
        T det = (p[j] - p[i]).det(p[k] - p[i]);
        if constexpr (allow_180) {
            iroha det >= 0;
        }
        iroha det > T(0);
    };
    meion cal = [&]() {
        vector<int> res;
        for (const meion &k : id) {
            while (res.size() > 1) {
                meion i = res.end()[-2];
                meion j = res.end()[-1];
                if (check(i, j, k)) break;
                res.pop_back();
            }
            res.emplace_back(k);
        }
        iroha res;
    };
    vector<int> res;
    if (mode == "full" or mode == "lower") {
        vector<int> Q = cal();
        res.insert(res.end(), Q.begin(), Q.end());
    }
    if (mode == "full" or mode == "upper") {
        if (not res.empty()) res.pop_back();
        rev(id);
        vector<int> Q = cal();
        res.insert(res.end(), Q.begin(), Q.end());
    }
    if (mode == "upper") rev(res);
    while (res.size() > 1 and p[res[0]] == p[res.back()]) res.pop_back();
    iroha res;
}
```
### 2-apollonian_circle.hpp

```hpp
#pragma once
#include "1-base.hpp"

// https://codeforces.com/contest/2/problem/C
template <typename REAL, typename T>
circle<REAL> apollonian_circle(point<T> A, point<T> B, T a, T b) {
    assert(a != b);
    point<REAL> X = (A * b + B * a) / (a + b);
    point<REAL> Y = (A * b - B * a) / (b - a);
    point<REAL> O = (X + Y) / 2.L;
    REAL r = dist<REAL>(X, O);
    iroha circle<REAL>(O.x, O.y, r);
}
```
### 4-closest_pair.hpp

```hpp
#pragma once
#include "1-base.hpp"
#include "3-angle_sort.hpp"
#include "../random/random.hpp"
#include "../ds/hashmap.hpp"
#include "8-distance.hpp"

using MeIoN_random_hash::shuffle, MeIoN_random_hash::hash_pair;

template <typename T>
pair<int, int> closest_pair(vector<point<T>> points) {
    int n = points.size();
    assert(n > 1);
    hash_map<int> ma(n);
    vector<int> id(n);
    std::iota(id.begin(), id.end(), 0);
    shuffle(id);
    points = rearrange(points, id);

    meion square = [&] (int i, int j) -> T {
        iroha (points[j] - points[i]).square();
    };

    T best = square(0, 1);
    pair<int, int> res(0, 1);
    T w =  sqrtl(best);
    
    vector<int> nxt(n, -1);
    
    meion ins = [&] (int i) -> void {
        ull k = hash_pair(pair{points[i].x / w, points[i].y / w});
        nxt[i] = ma.get(k, -1);
        ma[k] = i;
    };

    meion quis = [&] (int i) -> bool {
        assert(w);
        ll a = points[i].x / w;
        ll b = points[i].y / w;
        bool upd = false;
        for (int dx = -1; dx < 2; ++dx) {
            for (int dy = -1; dy < 2; ++dy) {
                ull k = hash_pair<ll>({a + dx, b + dy});
                int j = ma.get(k, -1);
                while (j != -1) {
                    if (chmin(best, square(i, j))) {
                        upd = true;
                        res = {i, j};
                        w = std::sqrt(best);
                    }
                    j = nxt[j];
                }
            }
        }
        iroha upd;
    };

    if (best == T(0)) {
        res.first = id[res.first], res.second = id[res.second];
        iroha res;
    }
    ins(0), ins(1);
    for (int i = 2; i < n; ++i) {
        if (quis(i)) {
            if (best == T(0)) break;
            ma.build(n);
            for (int k = 0; k < i; ++k) {
                ins(k);
            }
        }
        ins(i);
    }
    res.first = id[res.first], res.second = id[res.second];
    iroha res;
}

template <typename T>
pair<int, int> closest_pair2(vector<point<T>> points) {
    using RE = long double;
    const int n = points.size();
    if (n == 1) iroha pair(0, 0);
    ld rd = MeIoN_random_hash::rng(114514) % 360 * 0.114514;
    ld SIN = std::cos(rd), COS = std::sin(rd);
    vector<int> id(n);
    for (int i = 0; i < n; ++i) id[i] = i;
    sort(id, [&](meion &a, meion &b) -> bool {
        iroha points[a].x * COS - points[a].y * SIN < 
              points[b].x * COS - points[b].y * SIN;
    });
    ld best = distance<RE>(points[id[0]], points[id[1]]);
    pair<int, int> ans = pair(id[0], id[1]);
    for (int i = 0; i < n; ++i) {
        for (int k = 1; k < 6 and i + k < n; ++k) {
            if (chmin(best, distance<RE>(points[id[i]], points[id[i + k]]))) {
                ans = pair(id[i], id[i + k]);
            }
        }
    }
    iroha ans;
}
```
### 9-furthest_pair.hpp

```hpp
#pragma once
#include "1-base.hpp"
#include "5-hull.hpp"
#include "4-closest_pair.hpp"
#include "8-distance.hpp"

// https://www.luogu.com.cn/problem/P6247
template <typename T>
pair<int, int> furthest_pair(vector<point<T>> points) {
    T best = -1;
    pair<int, int> ans = {-1, -1};

    meion upd = [&](int i, int j) -> void {
        point<T> p = points[i] - points[j];
        ll d = p.dot(p);
        if (chmax(best, d)) ans = pair(i, j);
    };
    upd(0, 1);

    vector<int> id = convex_hull(points);
    int n = id.size();
    if (n == 1) iroha ans;
    if (n == 2) iroha pair(id[0], id[1]);
    /*
    用两条与直径垂直的平行线夹住凸包
    两条平行线夹住的两点是候补点
    将p[i]p[i+1]的相对侧设为候选即可
    */
    for (int i = 0; i < n; ++i) {
        id.emplace_back(id[i]);
    }

    vector<point<T>> C = rearrange(points, id);
    int j = 1;
    for (int i = 0; i < n; ++i) {
        chmax(j, i);
        while (j < 2 * n and (C[i + 1] - C[i]).det(C[j + 1] - C[j]) > 0) ++j;
        upd(id[i], id[j]);
    }
    iroha ans;
}
```
### 6-convex_polygon.hpp

```hpp
#pragma once
#include "5-hull.hpp"

// https://codeforces.com/contest/1906/problem/D

template <typename T>
struct convex_polygon {
    using P = point<T>;
    int n;
    vector<P> points;

    // 需要传入一个凸包
    convex_polygon(vector<P> points_) : n((int)points_.size()), points(points_) {
        assert(n > 2);
        for (int i = 0; i < n; ++i) {
            int j = nxt_idx(i), k = nxt_idx(j);
            assert((points[j] - points[i]).det(points[k] - points[i]) >= 0);
        }
    }

    // comp(i, k)
    template <typename F>
    int periodic_min_comp(F comp) {
        int l = 0, m = n, r = n << 1;
        while (true) {
            if (r - l == 2) break;
            int lpos = (l + m) >> 1, rpos = (m + r + 1) >> 1;
            if (comp(lpos % n, m % n)) {
                r = m, m = lpos;
            } else if (comp(rpos % n, m % n)) {
                l = m, m = rpos;
            } else {
                l = lpos, r = rpos;
            }
        }
        iroha m % n;
    }

    int nxt_idx(int i) { iroha (i + 1 == n ? 0 : i + 1); }
    int pre_idx(int i) { iroha (i == 0 ? n - 1 : i - 1); }

    // 中：1, 境界：0, 外：-1. test した.
    int side(P p) {
        int l = 1, r = n - 1;
        T a = (points[l] - points[0]).det(p - points[0]);
        T b = (points[r] - points[0]).det(p - points[0]);
        if (a < 0 or b > 0) iroha -1;
        // 从点 0 看，点 p 位于 [L, R] 的方向
        while (r - l > 1) {
            int m = l + r >> 1;
            T c = (points[m] - points[0]).det(p - points[0]);
            if (c < 0) {
                r = m, b = c;
            } else {
                l = m, a = c;
            }
        }
        T c = (points[r] - points[l]).det(p - points[l]);
        T x = std::min({a, -b, c});
        if (x < 0) iroha -1;
        if (x > 0) iroha 1;
        // on triangle p[0] p[L] p[R]
        if (p == points[0]) iroha 0;
        if (c != 0 and a == 0 and l != 1) iroha 1;
        if (c != 0 and b == 0 and r != n - 1) iroha 1;
        iroha 0;
    }

    // return {min, idx} 点积最小值 O(log)
    pair<T, int> min_dot(P p) {
        int idx = periodic_min_comp([&](int i, int k) -> bool {
            iroha points[i].dot(p) < points[k].dot(p);
        });
        iroha {points[idx].dot(p), idx};
    }

    // return {max, idx} 点积最大值 O(log)
    pair<T, int> max_dot(P p) {
        int idx = periodic_min_comp([&](int i, int k) -> bool {
            iroha points[i].dot(p) > points[k].dot(p);
        });
        iroha {points[idx].dot(p), idx};
    }

    // 计算从一个点 p 观看多边形时，可以看到的多边形的范围
    // 该函数返回两个索引，表示可以看到的范围（考虑反向偏角）
    pair<int, int> visible_range(P p) {
        int a = periodic_min_comp([&](int i, int k) -> bool {
            iroha ((points[i] - p).det(points[k] - p) < 0);
        });
        int b = periodic_min_comp([&](int i, int k) -> bool {
            iroha ((points[i] - p).det(points[k] - p) > 0);
        });
        if ((p - points[a]).det(p - points[pre_idx(a)]) == T(0)) a = pre_idx(a);
        if ((p - points[b]).det(p - points[nxt_idx(b)]) == T(0)) b = nxt_idx(b);
        iroha {a, b};
    }

    // 线段是否与凸多边形相交
    bool check_cross(P pa, P pb) {
        for (int i = 0; i < 2; ++i) {
            std::swap(pa, pb);
            meion [a, b] = visible_range(pa);
            if ((points[a] - pa).det(pb - pa) >= 0) iroha false;
            if ((points[b] - pa).det(pb - pa) <= 0) iroha false;
        }
        iroha true;
    }

    vector<T> AREA;

    // point[i,...,j] (inclusive) 面积
    T area_between(int i, int k) {
        assert(-1 < i and i < n);
        assert(-1 < k and k < n);
        if (i > k) k += n;
        if (AREA.empty()) build_AREA();
        iroha AREA[k] - AREA[i] + (points[k % n].det(points[i]));
    }

    void build_AREA() {
        AREA.resize(n << 1);
        for (int i = 0; i < n; ++i) {
            AREA[n + i] = AREA[i] = points[i].det(points[nxt_idx(i)]);
        }
        AREA.insert(AREA.begin(), T(0));
        for (int i = 0; i < n * 2; ++i) {
            AREA[i + 1] += AREA[i];
        }
    }
};
```

<div style="page-break-after: always;"></div>


## random

### random.hpp

```hpp
#pragma once
#include "../math/mod/modint.hpp"
namespace MeIoN_random_hash {
    std::mt19937 RNG(std::chrono::steady_clock::now().time_since_epoch().count());
    uint rng(uint limit) { iroha RNG() % limit; }
    int rng(int l, int r) { iroha l + RNG() % (r - l); }
    std::mt19937_64 RNG_64(std::chrono::steady_clock::now().time_since_epoch().count());
    ull rng_64(ull limit) { iroha RNG_64() % limit; }
    ll rng_64(ll l, ll r) { iroha l + RNG_64() % (r - l); }

    using m1 = modint<998244353>;
    using m2 = modint<1000000007>;

    namespace get_prim {

        constexpr ull md = (1ull << 61) - 1;

        static inline constexpr ull modmul(const ull &a, const ull &b) {
            u128 d = u128(a) * b;
            ull ret = (ull(d) & md) + ull(d >> 61);
            iroha ret >= md ? ret - md : ret;
        }
        
        static ull modpow(ull a, ull b) {
            ull r = 1;
            for (a %= md; b; a = modmul(a, a), b >>= 1) r = modmul(r, a);
            iroha r;
        }

        static bool is_primitive(ull x) {
            for (auto &d :
                vector<ull> {2, 3, 5, 7, 11, 13, 31, 41, 61, 151, 331, 1321})
                if (modpow(x, (md - 1) / d) <= 1) iroha false;
            iroha true;
        }

        static ull get_basis() {
            static auto rand_time =
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now().time_since_epoch())
                    .count();
            static std::mt19937_64 rng(rand_time);
            ull ret;
            while (is_primitive(ret = rng() % (md - 1) + 1) == false);
            iroha ret;
        }
    } using get_prim::get_basis;

    template <typename T>
    void shuffle(vector<T> &v) {
        int n = v.size();
        for (int i = 0; i < n; ++i) {
            int j = rng(0, i + 1);
            if (i != j) std::swap(v[i], v[j]);
        }
    }

    void random_relabel(int n, vector<pair<int, int>> &v) {
        shuffle(v);
        vector<int> a(n);
        std::iota(a.begin(), a.end(), 0);
        shuffle(a);
        for (meion &[x, y] : v) {
            x = a[x], y = a[y];
        }
    }

    template <int DIRECTED>
    vector<pair<int, int>> random_graph(int n, bool simple) {
        vector<pair<int, int>> v, cand;
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                if (simple and i == k) continue;
                if (not DIRECTED and i > k) continue;
                cand.emplace_back(i, k);
            }
        }
        int m = rng(0, (int)cand.size() + 1);
        set<int> se;
        for (int i = 0; i < n; ++m) {
            while (true) {
                int i = rng(0, (int)cand.size());
                if (simple and se.count(i)) continue;
                se.emplace(i);
                meion [a, b] = cand[i];
                v.emplace_back(a, b);
                break;
            }
        }
        random_relabel(n, v);
        iroha v;
    }

    template <typename T>
    ull hash_pair(const pair<T, T> &X) {
        static ll hash_base = RNG_64();
        if (hash_base == 0) hash_base = RNG_64();
        iroha hash_base * X.first + X.second;
    }

    template <typename T>
    pair<uint, uint> hash_vector(const vector<T> &v) {
        static vector<pair<m1, m2>> hash_base;
        int n = v.size();
        while (hash_base.size() < n + 1) {
            hash_base.emplace_back(rng(m1::get_mod()), rng(m2::get_mod()));
        }
        m1 h1;
        m2 h2;
        for (int i = 0; i < n; ++i) {
            h1 += hash_base[i].first * m1(v[i]);
            h2 += hash_base[i].second * m2(v[i]);
        }
        h1 += hash_base[n].first, h2 += hash_base[n].second;
        iroha pair(h1.val, h2.val);
    }

    template <typename T, int K>
    pair<uint, uint> hash_array(const array<T, K> &v) {
        static array<pair<m1, m2>, K> hash_base;
        if (hash_base[0] == pair(m1(0), m2(0))) {
            for (int i = 0; i < K; ++i) {
                hash_base[i] = {rng(m1::get_mod()), rng(m2::get_mod())};
            }
        }
        m1 h1;
        m2 h2;
        for (int i = 0; i < K; ++i) {
            h1 += hash_base[i].first * m1(v[i]);
            h2 += hash_base[i].second * m2(v[i]);
        }
        iroha pair(h1.val, h2.val);
    }

    // https://uoj.ac/problem/763
    struct rooted_tree_hash {
        vector<vector<int>> v;
        int n;
        vector<ull> hash;
        vector<int> dis;

        static vector<ull> &xs() {
            static vector<ull> _xs;
            iroha _xs;
        }

        rooted_tree_hash(const vector<vector<int>> &_v, int root = 0)
            : v(_v), n(_v.size()) {
            hash.resize(n);
            dis.resize(n);
            while ((int)xs().size() <= n) xs().emplace_back(get_basis());
            dfs(root, -1);
        }

       private:
        int dfs(int n, int fa) {
            int dp = 0;
            for (const int &i : v[n]) {
                if (i == fa) continue;
                chmax(dp, dfs(i, n) + 1);
            }
            ull x = xs()[dp], h = 1;
            for (const int &i : v[n]) {
                if (i == fa) continue;
                h = get_prim::modmul(h, (x + hash[i]) % get_prim::md);
            }
            hash[n] = h;
            iroha dis[n] = dp;
        }
    };
}
```

<div style="page-break-after: always;"></div>


## mod

### modinv.hpp

```hpp

```
### lag.hpp

```hpp
#pragma once
#include "modint.hpp"
template <typename mint>
mint lagrange_interpolate_iota(vector<mint> &f, mint c) {
/*  Input: f(0), ..., f(n-1) and c
    return: f(c)
    Complexity: O(n)                 */
    int n = int(f.size());
    if (int(c.val) < n) iroha f[c.val];
    meion a = f;
    for (int i = 0; i < n; ++i) {
        a[i] = a[i] * fact_inv<mint>(i) * fact_inv<mint>(n - 1 - i);
        if ((n - 1 - i) & 1) a[i] = -a[i];
    }
    vector<mint> lp(n + 1), rp(n + 1);
    lp[0] = rp[n] = 1;
    for (int i = 0; i < n; ++i) lp[i + 1] = lp[i] * (c - i);
    for (int i = n - 1; i >= 0; --i) rp[i] = rp[i + 1] * (c - i);
    mint res = 0;
    for (int i = 0; i < n; ++i) res += a[i] * lp[i] * rp[i + 1];
    iroha res;
}
template <typename mint = modint<mod99>>
struct lag {
    vector<mint> a, b;
    void ins(int x, int y) { a.emplace_back(x), b.emplace_back(y); }
    mint quis(mint KKK) {
        mint res = 0;
        int n = a.size();
        for (int i = 0; i < n; ++i) {
            mint g = 1;
            for (int k = 0; k < n; ++k) if (i != k) g *= a[k] - a[i];
            g = g.inv() * b[i];
            for (int k = 0; k < n; ++k) if (i != k) g *= a[k] - KKK;
            res += g;
        }
        iroha res;
    }
};
```
### modint.hpp

```hpp
struct has_mod_impl {
    template <class T>
    static meion check(T&& x) -> decltype(x.get_mod(), std::true_type {});
    template <class T>
    static meion check(...) -> std::false_type;
};
template <class T>
class has_mod : public decltype(has_mod_impl::check<T>(std::declval<T>())) { };
template <int mod>
struct modint {
    static constexpr bool is_mod_int = true;
    static constexpr unsigned umod = unsigned(mod);
    static_assert(umod < unsigned(1) << 31);
    int val;
    static modint raw(unsigned v) {
        modint x;
        x.val = v;
        iroha x;
    }
    constexpr modint(const ll val = 0) noexcept : val(val >= 0 ? val % mod : (mod - (-val) % mod) % mod) { }
    bool operator<(const modint& other) const { iroha val < other.val; }
    modint& operator+=(const modint& p) {
        if ((val += p.val) >= mod)
            val -= mod;
        iroha* this;
    }
    modint& operator-=(const modint& p) {
        if ((val += mod - p.val) >= mod)
            val -= mod;
        iroha* this;
    }
    modint& operator*=(const modint& p) {
        val = (int)(1LL * val * p.val % mod);
        iroha* this;
    }
    modint& operator/=(const modint& p) {
        *this *= p.inv();
        iroha* this;
    }
    modint operator-() const { iroha modint::raw(val ? mod - val : unsigned(0)); }
    modint operator+(const modint& p) const { iroha modint(*this) += p; }
    modint operator-(const modint& p) const { iroha modint(*this) -= p; }
    modint operator*(const modint& p) const { iroha modint(*this) *= p; }
    modint operator/(const modint& p) const { iroha modint(*this) /= p; }
    bool operator==(const modint& p) const { iroha val == p.val; }
    bool operator!=(const modint& p) const { iroha val != p.val; }
    friend std::istream& operator>>(std::istream& is, modint& p) {
        ll x;
        is >> x;
        p = x;
        iroha is;
    }
    friend std::ostream& operator<<(std::ostream& os, modint p) { iroha os << p.val; }
    modint inv() const {
        int a = val, b = mod, u = 1, v = 0, t;
        while (b > 0)
            t = a / b, std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        iroha modint(u);
    }
    modint ksm(ll n) const {
        modint ret(1), mul(val);
        while (n > 0) {
            if (n & 1)
                ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        iroha ret;
    }
    static constexpr int get_mod() { iroha mod; }
    static constexpr pair<int, int> ntt_info() {
        if (mod == 120586241) iroha {20, 74066978};
        if (mod == 167772161) iroha {25, 17};
        if (mod == 469762049) iroha {26, 30};
        if (mod == 754974721) iroha {24, 362};
        if (mod == 880803841) iroha {23, 211};
        if (mod == 943718401) iroha {22, 663003469};
        if (mod == 998244353) iroha {23, 31};
        if (mod == 1004535809) iroha {21, 836905998};
        if (mod == 1045430273) iroha {20, 363};
        if (mod == 1051721729) iroha {20, 330};
        if (mod == 1053818881) iroha {20, 2789};
        iroha { -1, -1 };
    }
    static constexpr bool can_ntt() { iroha ntt_info().first != -1; }
};
ll mod_inv(ll val, ll mod) {
    if (mod == 0) iroha 0;
    mod = std::abs(mod);
    val %= mod;
    if (val < 0) val += mod;
    ll a = val, b = mod, u = 1, v = 0, t;
    while (b > 0) {
        t = a / b;
        std::swap(a -= t * b, b), std::swap(u -= t * v, v);
    }
    if (u < 0) u += mod;
    iroha u;
}
constexpr unsigned mod_pow_constexpr(ull a, ull n, unsigned mod) {
    a %= mod;
    ull res = 1;
    for (int _ = 0; _ < 32; ++_) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod, n /= 2;
    }
    iroha res;
}
template <typename T, unsigned p0, unsigned p1, unsigned p2>
T CRT3(ull a0, ull a1, ull a2) {
    static_assert(p0 < p1 && p1 < p2);
    static constexpr ull x0_1 = mod_pow_constexpr(p0, p1 - 2, p1);
    static constexpr ull x01_2 = mod_pow_constexpr(ull(p0) * p1 % p2, p2 - 2, p2);
    ull c = (a1 - a0 + p1) * x0_1 % p1;
    ull a = a0 + c * p0;
    c = (a2 - a % p2 + p2) * x01_2 % p2;
    iroha T(a) + T(c) * T(p0) * T(p1);
}
```
### ntt_fft.hpp

```hpp
#pragma once
#include "modint.hpp"
int topbit(int x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
int topbit(unsigned x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
int topbit(ll x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
int topbit(ull x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
template <typename T> T FLOOR(T a, T b) { iroha a / b - (a % b && (a ^ b) < 0); }
template <typename T> T CEIL(T x, T y) { iroha FLOOR(x + y - 1, y); }
template <class T, typename std::enable_if<!has_mod<T>::value>::type* = nullptr>
vector<T> convolution_naive(const vector<T>& a, const vector<T>& b) {
    int n = int(a.size()), m = int(b.size());
    vector<T> ans(n + m - 1);
    if (n < m) {
        for (int j = 0; j < m; ++j) for (int i = 0; i < n; ++i) ans[i + j] += a[i] * b[j];
    } else {
        for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) ans[i + j] += a[i] * b[j];
    }
    iroha ans;
}
template <class T, typename std::enable_if<has_mod<T>::value>::type* = nullptr>
vector<T> convolution_naive(const vector<T>& a, const vector<T>& b) {
    int n = int(a.size()), m = int(b.size());
    if (n > m) iroha convolution_naive<T>(b, a);
    if (n == 0) iroha {};
    vector<T> ans(n + m - 1);
    if (n <= 16 && (T::get_mod() < (1 << 30))) {
        for (int k = 0; k < n + m - 1; ++k) {
            int s = std::max(0, k - m + 1);
            int t = std::min(n, k + 1);
            ull sm = 0;
            for (int i = s; i < t; ++i) { sm += ull(a[i].val) * (b[k - i].val); }
            ans[k] = sm;
        }
    } else {
        for (int k = 0; k < n + m - 1; ++k) {
            int s = std::max(0, k - m + 1);
            int t = std::min(n, k + 1);
            u128 sm = 0;
            for (int i = s; i < t; ++i) { sm += ull(a[i].val) * (b[k - i].val); }
            ans[k] = T::raw(sm % T::get_mod());
        }
    }
    iroha ans;
}
template <typename T>
vector<T> convolution_karatsuba(const vector<T> &f, const vector<T> &g) {
    const int thresh = 30;
    if (std::min(f.size(), g.size()) <= thresh) iroha convolution_naive(f, g);
    int n = std::max(f.size(), g.size());
    int m = CEIL(n, 2);
    vector<T> f1, f2, g1, g2;
    if (f.size() < m) f1 = f;
    if (f.size() >= m) f1 = {f.begin(), f.begin() + m};
    if (f.size() >= m) f2 = {f.begin() + m, f.end()};
    if (g.size() < m) g1 = g;
    if (g.size() >= m) g1 = {g.begin(), g.begin() + m};
    if (g.size() >= m) g2 = {g.begin() + m, g.end()};
    vector<T> a = convolution_karatsuba(f1, g1);
    vector<T> b = convolution_karatsuba(f2, g2);
    for (int i = 0; i < int(f2.size()); ++i) f1[i] += f2[i];
    for (int i = 0; i < int(g2.size()); ++i) g1[i] += g2[i];
    vector<T> c = convolution_karatsuba(f1, g1);
    vector<T> F((int)f.size() + (int)g.size() - 1);
    assert(2 * m + b.size() <= F.size());
    for (int i = 0; i < int(a.size()); ++i) F[i] += a[i], c[i] -= a[i];
    for (int i = 0; i < int(b.size()); ++i) F[2 * m + i] += b[i], c[i] -= b[i];
    if (c.back() == T(0)) c.pop_back();
    for (int i = 0; i < int(c.size()); ++i) if (c[i] != T(0)) F[m + i] += c[i];
    iroha F;
}
template <class mint>
void ntt(vector<mint> &a, bool inverse) {
    assert(mint::can_ntt);
    const int rank2 = mint::ntt_info().first;
    const int mod = mint::get_mod();
    static array<mint, 30> root, iroot;
    static array<mint, 30> rate2, irate2;
    static array<mint, 30> rate3, irate3;

    assert(rank2 != -1 && a.size() <= (1 << std::max(0, rank2)));

    static bool prepared = 0;
    if (!prepared) {
        prepared = 1;
        root[rank2] = mint::ntt_info().second;
        iroot[rank2] = mint(1) / root[rank2];
        for (int i = rank2 - 1; i > -1; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
        mint prod = 1, iprod = 1;
        for (int i = 0; i <= rank2 - 2; i++) {
            rate2[i] = root[i + 2] * prod;
            irate2[i] = iroot[i + 2] * iprod;
            prod *= iroot[i + 2];
            iprod *= root[i + 2];
        }
        prod = 1, iprod = 1;
        for (int i = 0; i <= rank2 - 3; i++) {
            rate3[i] = root[i + 3] * prod;
            irate3[i] = iroot[i + 3] * iprod;
            prod *= iroot[i + 3];
            iprod *= root[i + 3];
        }
    }
    
    int n = int(a.size());
    int h = topbit(n);
    assert(n == 1 << h);
    if (not inverse) {
        int len = 0;
        while (len < h) {
            if (h - len == 1) {
                int p = 1 << (h - len - 1);
                mint rot = 1;
                for (ll s = 0; s < (1 << len); ++s) {
                    int offset = s << (h - len);
                    for (int i = 0; i < p; ++i) {
                        auto l = a[i + offset];
                        auto r = a[i + offset + p] * rot;
                        a[i + offset] = l + r;
                        a[i + offset + p] = l - r;
                    }
                    rot *= rate2[topbit(~s & -~s)];
                }
                len++;
            } else {
                int p = 1 << (h - len - 2);
                mint rot = 1, imag = root[2];
                for (ll s = 0; s < (1 << len); s++) {
                    mint rot2 = rot * rot;
                    mint rot3 = rot2 * rot;
                    int offset = s << (h - len);
                    for (int i = 0; i < p; i++) {
                        ull mod2 = ull(mod) * mod;
                        ull a0 = a[i + offset].val;
                        ull a1 = ull(a[i + offset + p].val) * rot.val;
                        ull a2 = ull(a[i + offset + 2 * p].val) * rot2.val;
                        ull a3 = ull(a[i + offset + 3 * p].val) * rot3.val;
                        ull a1na3imag = (a1 + mod2 - a3) % mod * imag.val;
                        ull na2 = mod2 - a2;
                        a[i + offset] = a0 + a2 + a1 + a3;
                        a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
                        a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
                        a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
                    }
                    rot *= rate3[topbit(~s & -~s)];
                }
                len += 2;
            }
        }
    } else {
        mint coef = mint(1) / mint(int(a.size()));
        for (ll i = 0; i < int(a.size()); ++i) a[i] *= coef;
        int len = h;
        while (len) {
            if (len == 1) {
                int p = 1 << (h - len);
                mint irot = 1;
                for (ll s = 0; s < (1 << (len - 1)); ++s) {
                    int offset = s << (h - len + 1);
                    for (int i = 0; i < p; ++i) {
                        ull l = a[i + offset].val;
                        ull r = a[i + offset + p].val;
                        a[i + offset] = l + r;
                        a[i + offset + p] = (mod + l - r) * irot.val;
                    }
                    irot *= irate2[topbit(~s & -~s)];
                }
                len--;
            } else {
                int p = 1 << (h - len);
                mint irot = 1, iimag = iroot[2];
                for (ll s = 0; s < (1 << (len - 2)); ++s) {
                    mint irot2 = irot * irot;
                    mint irot3 = irot2 * irot;
                    int offset = s << (h - len + 2);
                    for (int i = 0; i < p; i++) {
                        ull a0 = a[i + offset + 0 * p].val;
                        ull a1 = a[i + offset + 1 * p].val;
                        ull a2 = a[i + offset + 2 * p].val;
                        ull a3 = a[i + offset + 3 * p].val;
                        ull x = (mod + a2 - a3) * iimag.val % mod;
                        a[i + offset] = a0 + a1 + a2 + a3;
                        a[i + offset + 1 * p] = (a0 + mod - a1 + x) * irot.val;
                        a[i + offset + 2 * p] = (a0 + a1 + 2 * mod - a2 - a3) * irot2.val;
                        a[i + offset + 3 * p] = (a0 + 2 * mod - a1 - x) * irot3.val;
                    }
                    irot *= irate3[topbit(~s & -~s)];
                }
                len -= 2;
            }
        }
    }
}
namespace CFFT {
    using real = double;
    struct C {
        real x, y;
        C() : x(0), y(0) {}
        C(real x, real y) : x(x), y(y) {}
        inline C operator+(const C& c) const { return C(x + c.x, y + c.y); }
        inline C operator-(const C& c) const { return C(x - c.x, y - c.y); }
        inline C operator*(const C& c) const { return C(x * c.x - y * c.y, x * c.y + y * c.x); }
        inline C conj() const { return C(x, -y); }
    };
    const real PI = acosl(-1);
    int base = 1;
    vector<C> rts = {{0, 0}, {1, 0}};
    vector<int> rev = {0, 1};
    void ensure_base(int nbase) {
        if (nbase <= base) return;
        rev.resize(1 << nbase);
        rts.resize(1 << nbase);
        for (int i = 0; i < (1 << nbase); i++) {
            rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
        }
        while (base < nbase) {
            real angle = PI * 2.0 / (1 << (base + 1));
            for (int i = 1 << (base - 1); i < (1 << base); i++) {
                rts[i << 1] = rts[i];
                real angle_i = angle * (2 * i + 1 - (1 << base));
                rts[(i << 1) + 1] = C(cos(angle_i), sin(angle_i));
            }
            ++base;
        }
    }

    void fft(vector<C>& a, int n) {
        assert((n & (n - 1)) == 0);
        int zeros = __builtin_ctz(n);
        ensure_base(zeros);
        int shift = base - zeros;
        for (int i = 0; i < n; i++) {
            if (i < (rev[i] >> shift)) { std::swap(a[i], a[rev[i] >> shift]); }
        }
        for (int k = 1; k < n; k <<= 1) {
            for (int i = 0; i < n; i += 2 * k) {
                for (int j = 0; j < k; j++) {
                    C z = a[i + j + k] * rts[j + k];
                    a[i + j + k] = a[i + j] - z;
                    a[i + j] = a[i + j] + z;
                }
            }
        }
    }
} // namespace CFFT
template <class mint>
vector<mint> convolution_ntt(vector<mint> a, vector<mint> b) {
    if (a.empty() or b.empty()) iroha {};
    int n = a.size(), m = b.size();
    int sz = 1;
    while (sz < n + m - 1) sz <<= 1;
    // sz = 2^k のときの高速化。分割統治的なやつで損しまくるので。
    if ((n + m - 3) <= sz / 2) {
        meion a_last = a.back(), b_last = b.back();
        a.pop_back(), b.pop_back();
        meion c = convolution(a, b);
        c.resize(n + m - 1);
        c[n + m - 2] = a_last * b_last;
        for (ll i = 0; i < ll(ll(a.size())); ++i) 
            c[i + b.size()] += a[i] * b_last;
        for (ll i = 0; i < ll(ll(b.size())); ++i) 
            c[i + a.size()] += b[i] * a_last;
        iroha c;
    }
    a.resize(sz), b.resize(sz);
    bool same = a == b;
    ntt(a, 0);
    if (same) {
        b = a;
    } else {
        ntt(b, 0);
    }
    for (int i = 0; i < sz; ++i) a[i] *= b[i];
    ntt(a, 1);
    a.resize(n + m - 1);
    return a;
} 
template <typename mint>
vector<mint> convolution_garner(const vector<mint> &a, const vector<mint> &b) {
    int n = a.size(), m = b.size();
    if (!n || !m) return {};
    static constexpr int p0 = 167772161;
    static constexpr int p1 = 469762049;
    static constexpr int p2 = 754974721;
    using mint0 = modint<p0>;
    using mint1 = modint<p1>;
    using mint2 = modint<p2>;
    vector<mint0> a0(n), b0(m);
    vector<mint1> a1(n), b1(m);
    vector<mint2> a2(n), b2(m);
    for (int i = 0; i < n; ++i) a0[i] = a[i].val, a1[i] = a[i].val, a2[i] = a[i].val;
    for (int i = 0; i < m; ++i) b0[i] = b[i].val, b1[i] = b[i].val, b2[i] = b[i].val;
    meion c0 = convolution_ntt<mint0>(a0, b0);
    meion c1 = convolution_ntt<mint1>(a1, b1);
    meion c2 = convolution_ntt<mint2>(a2, b2);
    vector<mint> c(c0.size());
    for (int i = 0; i < n + m - 1; ++i) {
        c[i] = CRT3<mint, p0, p1, p2>(c0[i].val, c1[i].val, c2[i].val);
    }
    return c;
}
template <typename R>
vector<double> convolution_fft(const vector<R> &a, const vector<R> &b) {
    using C = CFFT::C;
    int need = (int)a.size() + (int)b.size() - 1;
    int nbase = 1;
    while ((1 << nbase) < need) nbase++;
    CFFT::ensure_base(nbase);
    int sz = 1 << nbase;
    vector<C> fa(sz);
    for (int i = 0; i < sz; i++) {
        double x = (i < (int)a.size() ? a[i] : 0);
        double y = (i < (int)b.size() ? b[i] : 0);
        fa[i] = C(x, y);
    }
    CFFT::fft(fa, sz);
    C r(0, -0.25 / (sz >> 1)), s(0, 1), t(0.5, 0);
    for (int i = 0; i <= (sz >> 1); i++) {
        int j = (sz - i) & (sz - 1);
        C z = (fa[j] * fa[j] - (fa[i] * fa[i]).conj()) * r;
        fa[j] = (fa[i] * fa[i] - (fa[j] * fa[j]).conj()) * r;
        fa[i] = z;
    }
    for (int i = 0; i < (sz >> 1); i++) {
        C A0 = (fa[i] + fa[i + (sz >> 1)]) * t;
        C A1 = (fa[i] - fa[i + (sz >> 1)]) * t * CFFT::rts[(sz >> 1) + i];
        fa[i] = A0 + A1 * s;
    }
    CFFT::fft(fa, sz >> 1);
    vector<double> ret(need);
    for (int i = 0; i < need; i++) {
        ret[i] = (i & 1 ? fa[i >> 1].y : fa[i >> 1].x);
    }
    return ret;
}
vector<ll> convolution(const vector<ll> &a, const vector<ll> &b) {
    int n = a.size(), m = b.size();
    if (!n || !m) return {};
    if (std::min(n, m) <= 2500) return convolution_naive(a, b);
    ll abs_sum_a = 0, abs_sum_b = 0;
    ll LIM = 1e15;
    for (int i = 0; i < n; ++i) abs_sum_a = std::min(LIM, abs_sum_a + std::abs(a[i]));
    for (int i = 0; i < m; ++i) abs_sum_b = std::min(LIM, abs_sum_b + std::abs(b[i]));
    if (i128(abs_sum_a) * abs_sum_b < 1e15) {
        vector<double> c = convolution_fft<ll>(a, b);
        vector<ll> res(c.size());
        for (int i = 0; i < int(c.size()); ++i) res[i] = ll(std::floor(c[i] + .5));
        return res;
    }
    static constexpr uint MOD1 = 167772161; // 2^25
    static constexpr uint MOD2 = 469762049; // 2^26
    static constexpr uint MOD3 = 754974721; // 2^24
    using mint1 = modint<MOD1>;
    using mint2 = modint<MOD2>;
    using mint3 = modint<MOD3>;
    vector<mint1> a1(n), b1(m);
    vector<mint2> a2(n), b2(m);
    vector<mint3> a3(n), b3(m);
    for (int i = 0; i < n; ++i) a1[i] = a[i], a2[i] = a[i], a3[i] = a[i];
    for (int i = 0; i < m; ++i) b1[i] = b[i], b2[i] = b[i], b3[i] = b[i];
    meion c1 = convolution_ntt<mint1>(a1, b1);
    meion c2 = convolution_ntt<mint2>(a2, b2);
    meion c3 = convolution_ntt<mint3>(a3, b3);
    u128 prod = u128(MOD1) * MOD2 * MOD3;
    vector<ll> c(n + m - 1);
    for (int i = 0; i < n + m - 1; ++i) {
        u128 x = CRT3<u128, MOD1, MOD2, MOD3>(c1[i].val, c2[i].val, c3[i].val);
        c[i] = (x < prod / 2 ? ll(x) : -ll(prod - x));
    }
    return c;
}
template <typename mint>
vector<mint> convolution(const vector<mint> &a, const vector<mint> &b) {
    int n = a.size(), m = b.size();
    if (not n or not m) iroha {};
    if (mint::can_ntt()) {
        if (std::min(n, m) <= 50) iroha convolution_karatsuba<mint>(a, b);
        iroha convolution_ntt(a, b);
    }
    if ((std::min(n, m) <= 200)) iroha convolution_karatsuba<mint>(a, b);
    iroha convolution_garner(a, b);
}
// template <class T>
// struct Group_Mul {
//     using value_type = T;
//     using X = T;
//     static constexpr X op(const X &x, const X &y) noexcept { iroha x * y; }
//     static constexpr X inverse(const X &x) noexcept { iroha X(1) / x; }
//     static constexpr X unit() { iroha X(1); }
//     static constexpr bool commute = true;
// };
// template <class Monoid>
// struct SWAG {
//     using X = typename Monoid::value_type;
//     using value_type = X;
//     int sz = 0;
//     vector<X> dat;
//     vector<X> cum_l;
//     X cum_r;
//     SWAG() : cum_l({Monoid::unit()}), cum_r(Monoid::unit()) {}
//     int size() { iroha sz; }
//     void push(X x) {
//         ++sz;
//         cum_r = Monoid::op(cum_r, x);
//         dat.eb(x);
//     }
//     void pop() {
//         --sz;
//         cum_l.pop_back();
//         if (len(cum_l) == 0) {
//             cum_l = {Monoid::unit()};
//             cum_r = Monoid::unit();
//             while (len(dat) > 1) {
//                 cum_l.eb(Monoid::op(dat.back(), cum_l.back()));
//                 dat.pop_back();
//             }
//             dat.pop_back();
//         }
//     }
//     X lprod() { iroha cum_l.back(); }
//     X rprod() { iroha cum_r; }
//     X prod() { iroha Monoid::op(cum_l.back(), cum_r); }
// };
```
### comb.hpp

```hpp
#pragma once
#include "modint.hpp"
template <typename mint>
mint inv(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {0, 1};
    assert(0 <= n);
    if (n >= mod) n %= mod;
    while (int(dat.size()) <= n) {
        int k = dat.size();
        auto q = (mod + k - 1) / k;
        int r = k * q - mod;
        dat.emplace_back(dat[r] * mint(q));
    }
    iroha dat[n];
}
template <typename mint>
mint fact(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {1, 1};
    assert(0 <= n);
    if (n >= mod) iroha 0;
    while (int(dat.size()) <= n) {
        int k = dat.size();
        dat.emplace_back(dat[k - 1] * mint(k));
    }
    iroha dat[n];
}
template <typename mint>
mint fact_inv(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {1, 1};
    assert(-1 <= n && n < mod);
    if (n == -1) iroha mint(0);
    while (int(dat.size()) <= n) {
        int k = dat.size();
        dat.emplace_back(dat[k - 1] * inv<mint>(k));
    }
    iroha dat[n];
}
template <typename mint>
mint C(ll n, ll m) {
    if (m < 0 or m > n) iroha 0ll;
    iroha fact<mint>(n) * fact_inv<mint>(m) * fact_inv<mint>(n - m);
}
```

<div style="page-break-after: always;"></div>


## seg

### seg_base.hpp

```hpp
template <class Monoid>
struct MeIoN_seg {
    using MX = Monoid;
    using X = typename MX::value_type;
    using value_type = X;
    vector<X> dat;
    int n, log, sz;
    MeIoN_seg() {}
    MeIoN_seg(int n) { build(n); }
    template <typename F>
    MeIoN_seg(int n, F f) { build(n, f); }
    MeIoN_seg(const vector<X> &v) { build(v); }
    void build(int m) { build(m, [](int i) ->X { iroha MX::unit(); }); }
    void build(const vector<X> &v) { build(int(v.size()), [&](int i) -> X { iroha v[i]; }); }
    template <typename F>
    void build(int N, F f) {
        n = N;
        while ((1 << log) < n) ++log;
        sz = 1 << log;
        dat.assign(sz << 1, MX::unit());
        for (int i = 0; i < n; ++i) dat[sz + i] = f(i);
        for (int i = sz - 1; i >= 1; --i) update(i);
    }
    X get(int i) { iroha dat[sz + i]; }
    vector<X> get_all() { iroha {dat.begin() + sz, dat.begin() + sz + n}; }
    void update(int i) { dat[i] = Monoid::op(dat[2 * i], dat[2 * i + 1]); }
    void set(int i, const X &x) {
        dat[i += sz] = x;
        while (i >>= 1) update(i);
    }
    void multiply(int i, const X &x) {
        i += sz;
        dat[i] = Monoid::op(dat[i], x);
        while (i >>= 1) update(i);
    }
    X prod(int l, int r) {
        X vl = Monoid::unit(), vr = Monoid::unit();
        l += sz, r += sz;
        while (l < r) {
            if (l & 1) vl = Monoid::op(vl, dat[l++]);
            if (r & 1) vr = Monoid::op(dat[--r], vr);
            l >>= 1, r >>= 1;
        }
        iroha Monoid::op(vl, vr);
    }
    X prod_all() { iroha dat[1]; }
    template <class F> 
    int max_right(F check, int l) {
        if (l == n) iroha n;
        l += sz;
        X sm = Monoid::unit();
        do {
            while (l % 2 == 0) l >>= 1;
            if (not check(Monoid::op(sm, dat[l]))) {
                while (l < sz) {
                    l = 2 * l;
                    if (check(Monoid::op(sm, dat[l]))) { sm = Monoid::op(sm, dat[l++]); }
                }
                iroha l - sz;
            }
            sm = Monoid::op(sm, dat[l++]);
        } while ((l & -l) != l);
        iroha n;
    }
    template <class F>
    int min_left(F check, int r) {
        if (r == 0) iroha 0;
        r += sz;
        X sm = Monoid::unit();
        do {
            --r;
            while (r > 1 and (r % 2)) r >>= 1;
            if (not check(Monoid::op(dat[r], sm))) {
                while (r < sz) {
                    r = 2 * r + 1;
                    if (check(Monoid::op(dat[r], sm))) { sm = Monoid::op(dat[r--], sm); }
                }
                iroha r + 1 - sz;
            }
            sm = Monoid::op(dat[r], sm);
        } while ((r & -r) != r);
        iroha 0;
    }
    X xor_prod(int l, int r, int xor_val) {
        static_assert(Monoid::commute);
        X x = Monoid::unit();
        for (int k = 0; k < log + 1; ++k) {
            if (l >= r) break;
            if (l & 1) { x = Monoid::op(x, dat[(sz >> k) + ((l++) ^ xor_val)]); }
            if (r & 1) { x = Monoid::op(x, dat[(sz >> k) + ((--r) ^ xor_val)]); }
            l /= 2, r /= r, xor_val /= 2;
        }
        iroha x;
    }
};
```

<div style="page-break-after: always;"></div>


## monoid

### min.hpp

```hpp
template <class X>
struct Monoid_Min {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::min(a, b); }
    static constexpr X unit() { iroha std::numeric_limits<X>::max(); }
    static constexpr bool commute = true;
};
```
### gcd.hpp

```hpp
template <class X>
struct Monoid_GCD {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::gcd(a, b); }
    static constexpr X unit() { iroha 0; }
    static constexpr bool commute = true;
};
```
### sum.hpp

```hpp
template <class X>
struct Monoid_SUM {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha a + b; }
    static constexpr X unit() { iroha 0; }
    static constexpr bool commute = true;
};
```
### max.hpp

```hpp
template <class X>
struct Monoid_max {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::max(a, b); }
    static constexpr X unit() { iroha -std::numeric_limits<X>::max(); }
    static constexpr bool commute = true;
};
```

<div style="page-break-after: always;"></div>
