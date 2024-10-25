# MeIoN's XCPC Template

<style>
h3 { page-break-before: avoid; }
</style>

## 目录

- [MeIoN's XCPC Template](#meions-xcpc-template)
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
    - [1-base.hpp](#1-basehpp)
    - [2-apollonian\_circle.hpp](#2-apollonian_circlehpp)
    - [3-angle\_sort.hpp](#3-angle_sorthpp)
    - [4-closest\_pair.hpp](#4-closest_pairhpp)
    - [5-hull.hpp](#5-hullhpp)
    - [6-convex\_polygon.hpp](#6-convex_polygonhpp)
    - [7-points\_in\_triangles.hpp](#7-points_in_triangleshpp)
    - [8-distance.hpp](#8-distancehpp)
    - [9-furthest\_pair.hpp](#9-furthest_pairhpp)
    - [10-triangle\_area.hpp](#10-triangle_areahpp)
    - [11-in\_circle.hpp](#11-in_circlehpp)
  - [random](#random)
    - [random.hpp](#randomhpp)
  - [mod](#mod)
    - [lag.hpp](#laghpp)
    - [modint.hpp](#modinthpp)




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
        const auto [n, fa] = st.back();
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
    ret.reserve(2);
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
    return ret;
}
```
### LTT.hpp

```hpp
vector<int> get_fa(const vector<vector<int>> &v, int s) {
    int n = v.size();
    vector<int> pos(n, -1), p, label(n), dom(n), sdom(n), dsu(n), par(n);
    vector<vector<int>> rg(n), bucket(n);
    auto dfs = [&] (auto &&se, int n)->void {
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
    auto find = [&] (auto &&se, int n, int x) {
        if (n == dsu[n]) return x ? -1 : n;
        int v = se(se, dsu[n], x + 1);
        if (v < 0) return n;
        if (sdom[label[dsu[n]]] < sdom[label[n]]) {
            label[n] = label[dsu[n]];
        }
        dsu[n] = v;
        return x ? v : label[n];
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
    return res;
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
    return res;
}
```
### LCA.hpp

```hpp
template <const int N> struct LCA {
public:
    LCA (vector<vector<int>> _v, int rt) : 
    sz(_v.size()), v(_v), root(rt), up(sz), dis(sz), lg(0) {
        for (auto &i : up) i.fill(0);
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
        if (x == y) return x;
        for (int i = lg - 1; ~i; i--) {
            int X = up[x][i], Y = up[y][i];
            if (X != Y) x = X, y = Y;
        }
        return up[x][0];
    }
    int dist(int x,int y){
        return dis[x] + dis[y] - 2 * dis[lca(x, y)];
    }
private:
    int root, sz, lg;
    std::vector<std::vector<int>> v;
    std::vector<std::array<int, N>> up;
    std::vector<int> dis;
    void dfs (int n, int fa, int dp) { dis[n] = dp; up[n][0] = fa; for(int i = 1; i <= lg - 1; i++) up[n][i] = up[up[n][i - 1]][i - 1]; for (const auto &x : v[n]) { if(x == fa) continue; dfs(x, n, dp + 1); } }
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
    static u64 reduce(u128 b) { return (b + u128(u64(b) * r) * m) >> 64; }
    u64 x;
    m64() : x(0) {}
    m64(u64 x) : x(reduce(u128(x) * n2)){};
    u64 val() const { u64 y = reduce(x); return y >= m ? y - m : y; }
    m64 &operator+=(m64 y) { x += y.x - (m << 1); x = (i64(x) < 0 ? x + (m << 1) : x); return *this; }
    m64 &operator-=(m64 y) { x -= y.x; x = (i64(x) < 0 ? x + (m << 1) : x); return *this; }
    m64 &operator*=(m64 y) { x = reduce(u128(x) * y.x); return *this; }
    m64 operator+(m64 y) const { return m64(*this) += y; }
    m64 operator-(m64 y) const { return m64(*this) -= y; }
    m64 operator*(m64 y) const { return m64(*this) *= y; }
    bool operator==(m64 y) const { return (x >= m ? x - m : x) == (y.x >= m ? y.x - m : y.x); }
    bool operator!=(m64 y) const { return not operator==(y); }
    m64 pow(u64 n) const { m64 y = 1, z = *this; for (; n; n >>= 1, z *= z) if (n & 1) y *= z; return y; }
};

bool primetest(const uint64_t x) {
    using u64 = uint64_t;
    if (x == 2 or x == 3 or x == 5 or x == 7) return true;
    if (x % 2 == 0 or x % 3 == 0 or x % 5 == 0 or x % 7 == 0) return false;
    if (x < 121) return x > 1;
    const u64 d = (x - 1) >> __builtin_ctzll(x - 1);
    m64::set_mod(x);
    const m64 one(1), minus_one(x - 1);
    auto ok = [&](u64 a) {
        auto y = m64(a).pow(d);
        u64 t = d;
        while (y != one and y != minus_one and t != x - 1) y *= y, t <<= 1;
        if (y != minus_one and t % 2 == 0) return false;
        return true;
    };
    if (x < (1ull << 32)) {
        for (u64 a: {2, 7, 61}) if (not ok(a)) return false;
    } else { 
        for (u64 a: {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) { if (x <= a) return true; if (not ok(a)) return false; } 
    }
    return true;
}
ll rho(ll n, ll c) {
    m64::set_mod(n);
    assert(n > 1);
    const m64 cc(c);
    auto f = [&](m64 x) { return x * x + cc; };
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
    return g;
}
std::mt19937_64 rng_mt{std::random_device{}()};
ll rnd(ll n) { return std::uniform_int_distribution<ll>(0, n - 1)(rng_mt); }
ll find_prime_factor(ll n) {
    assert(n > 1);
    if (primetest(n)) return n;
    for (ll _ = 0; _ < 100ll; ++_) {
        ll m = rho(n, rnd(n));
        if (primetest(m)) return m;
        n = m;
    }
    std::cerr << "failed" << std::endl;
    assert(false);
    return -1;
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
    return pf;
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
    return res;
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
        for (auto p : primes) {
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
        return *this = res; 
    }
    MT operator*(const MT& p) { 
        return MT(*this) *= p; 
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
            (*this) = res; return res;
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
        return res;
    }
};
```
### Big_int.hpp

```hpp
namespace BGI {
    class ZeroDivisionError : public std::exception {
    public:
        const char* what() const throw() {return "BigInteger::divmod";}
    };
    class FFTLimitExceededError : public std::exception {
    public:
        const char* what() const throw() {return "BigInteger::fft_mul";}
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
            return out;
        }
        friend std::istream& operator>> (std::istream& in, BigInteger& x) {
            std::string s; in >> s; x = s; 
            return in;
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
        if (flag && !x.flag) return 1;
        if (!flag && x.flag) return -1;
        
        int sgn = (flag && x.flag ? 1 : -1);
        if (size > x.size) return sgn;
        if (size < x.size) return -sgn;
        
        for (int i = size; i >= 1; i--) {
            if (digits[i] > x.digits[i]) return sgn;
            if (digits[i] < x.digits[i]) return -sgn;
        }
        return 0;
    }

    inline void BigInteger::reserve(const int& sz) {
        if (sz < 0) return;
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
        return *this;
    }
    BigInteger& BigInteger::operator= (const long long& x) {
        flag = (x >= 0), reserve(4);
        if (x == 0) return size = 1, digits[1] = 0, *this;
        if (x == (-9223372036854775807ll - 1)) return *this = "-9223372036854775808";
        long long n = std::abs(x);
        do {push(n % BASE), n /= BASE;} while (n);
        return *this;
    }
    BigInteger& BigInteger::operator= (const std::string& s) {
        flag = true, reserve(s.size() / WIDTH + 1);
        if (s.empty() || s == "-") return *this = 0;
        int i = 0; if (s[0] == '-') flag = false, i++;
        for (int j = s.size() - 1; j >= i; j -= WIDTH) {
            int start = std::max(i, j - WIDTH + 1), len = j - start + 1;
            push(std::stoll(s.substr(start, len)));
        }
        return *this;
    }

    BigInteger& BigInteger::operator= (const std::vector<bool>& b) {
        *this = 0;
        if (b.empty() || (b.size() == 1 && b[0] == 0)) return *this;
        BigInteger pow2 = 1;
        for (int i = b.size() - 1; i >= 0; i--, pow2 += pow2) if (b[i]) *this += pow2;
        return *this;
    }

    void BigInteger::clear() {if (digits != nullptr) delete[] digits, digits = nullptr;}

    std::string BigInteger::to_string() const {std::stringstream ss; ss << *this; return ss.str();}
    long long BigInteger::to_long_long() const {return std::stoll(to_string());}
    std::vector<bool> BigInteger::to_binary() const {
        if (*this == 0) return {0};
        std::vector<bool> res;
        for (BigInteger x = *this; x != 0; x = x.div2()) res.emplace_back(x.digits[1] & 1);
        std::reverse(res.begin(), res.end());
        return res;
    };

    BigInteger BigInteger::operator- () const {
        if (*this == 0) return 0;
        BigInteger res = *this; res.flag = !flag; return res;
    }
    BigInteger BigInteger::abs() const {BigInteger res = *this; res.flag = true; return res;}

    bool BigInteger::operator== (const BigInteger& x) const {return compare(x) == 0;}
    #if __cplusplus >= 202002L
    auto BigInteger::operator<=> (const BigInteger& x) const {return compare(x);}
    #else
    bool BigInteger::operator< (const BigInteger& x) const {return compare(x) < 0;}
    bool BigInteger::operator> (const BigInteger& x) const {return compare(x) > 0;}
    bool BigInteger::operator!= (const BigInteger& x) const {return compare(x) != 0;}
    bool BigInteger::operator<= (const BigInteger& x) const {return compare(x) <= 0;}
    bool BigInteger::operator>= (const BigInteger& x) const {return compare(x) >= 0;}
    #endif //__cplusplus >= 202002L

    BigInteger BigInteger::operator+ (const BigInteger& x) const {
        if (!x.flag) return *this - x.abs();
        if (!flag) return x - abs();
        
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
        return res;
    }
    BigInteger BigInteger::operator- (const BigInteger& x) const {
        if (!x.flag) return *this + x.abs();
        if (!flag) return -(abs() + x);
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
        return res;
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
            
            complex operator+ (const complex& other) const {return complex(real + other.real, imag + other.imag);}
            complex operator- (const complex& other) const {return complex(real - other.real, imag - other.imag);}
            complex operator* (const complex& other) const {return complex(real * other.real - imag * other.imag, real * other.imag + other.real * imag);}
            
            complex& operator+= (const complex& other) {return real += other.real, imag += other.imag, *this;}
            complex& operator-= (const complex& other) {return real -= other.real, imag -= other.imag, *this;}
            complex& operator*= (const complex& other) {return *this = *this * other;}
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
            if (n <= 1) return;
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
            if (n <= 1) return;
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
        return res;
    }

    BigInteger BigInteger::operator* (const BigInteger& x) const {
        BigInteger zero = 0;
        if (*this == zero || x == zero) return zero;
        int n = size, m = x.size;
        long long lim = 1LL * n * m;
        
        if (lim >= FFT_LIMIT) {
            BigInteger res = fft_mul(*this, x);
            res.flag = !(flag ^ x.flag);
            return res;
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
        return res;
    }

    BigInteger& BigInteger::operator*= (int x) {
        if (x == 0 || *this == 0) return *this = 0;
        if (x < 0) flag = !flag, x = -x;
        digit_t carry = 0;
        for (int i = 1; i <= size || carry; i++) {
            if (i > size) push(0);
            digit_t cur = digits[i] * x + carry;
            carry = cur / BigInteger::BASE;
            digits[i] = cur % BigInteger::BASE;
        }
        while (size > 1 && digits[size] == 0) pop();
        return *this;
    }
    BigInteger BigInteger::operator* (const int& x) const {BigInteger t = *this; return t *= x;}

    BigInteger BigInteger::div2() const {
        BigInteger res = *this;
        for (int i = size; i >= 1; i--) {
            if ((res.digits[i] & 1) && (i > 1)) res.digits[i - 1] += BASE;
            res.digits[i] >>= 1;
        }
        while (res.size > 1 && res.digits[res.size] == 0) res.pop();
        return res;
    }
    BigInteger BigInteger::operator/ (const long long& x) const {
        if (x == 0) throw -1;
        if (*this == 0) return 0;
        if (x == 2) return div2();
        if (x == -2) {BigInteger res = div2(); res.flag = !res.flag; return res;}
        
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
        return res;
    }

    inline BigInteger BigInteger::move_r(int d) const {
        if (*this == 0 || d >= size) return 0;
        if (d == 0) return *this;
        BigInteger res; res.reserve(size - d + 1);
        for (int i = d + 1; i <= size; i++) res.push(digits[i]);
        return res;
    }
    inline BigInteger BigInteger::move_l(int d) const {
        if (*this == 0) return 0;
        if (d == 0) return *this;
        BigInteger res; res.reserve(size + d + 1);
        for (int i = 1; i <= d; i++) res.push(0);
        for (int i = 1; i <= size; i++) res.push(digits[i]);
        return res;
    }

    BigInteger BigInteger::newton_inv(int n) const {
        if (*this == 0) throw ZeroDivisionError();
        if (std::min(size, n - size) <= NEWTON_MIN_LEVEL) {
            BigInteger a; a.resize(n + 1);
            std::memset(a.digits, 0, sizeof(digit_t) * a.size);
            a.digits[n + 1] = 1;
            return a.divmod(*this, true).first;
        }
        int k = (n - size + 2) >> 1, k2 = k > size ? 0 : size - k;
        BigInteger x = move_r(k2);
        int n2 = k + x.size;
        BigInteger y = x.newton_inv(n2), a = y + y, b = (*this) * y * y;
        
        return a.move_l(n - n2 - k2) - b.move_r(2 * (n2 + k2) - n) - 1;
    }

    std::pair<BigInteger, BigInteger> BigInteger::newton_div(const BigInteger& x) const {
        int k = size - x.size + 2, k2 = k > x.size ? 0 : x.size - k;
        BigInteger x2 = x.move_r(k2);
        if (k2 != 0) x2 += 1;
        int n2 = k + x2.size;
        BigInteger u = (*this) * x2.newton_inv(n2);
        BigInteger q = u.move_r(n2 + k2), r = (*this) - q * x;
        while (r >= x) q += 1, r -= x;
        return std::make_pair(q, r);
    }

    std::pair<BigInteger, BigInteger> BigInteger::divmod(const BigInteger& x, bool dis_newton) const {
        static const int base = BigInteger::BASE;
        BigInteger a = abs(), b = x.abs();
        if (b == 0) throw ZeroDivisionError();
        if (a < b) return std::make_pair(0, flag ? a : -a);
        if (!dis_newton && size > NEWTON_LIMIT) return newton_div(x);
        
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
        return std::make_pair(q, r / t);
    }
    BigInteger BigInteger::operator/ (const BigInteger& x) const {return divmod(x).first;}

    BigInteger BigInteger::operator% (const long long& x) const {
        if (x == 2) return digits[1] & 1;
        if (x == 5) return digits[1] % 5;
        return *this - (*this / x * x);
    } 
    BigInteger BigInteger::operator% (const BigInteger& x) const {return divmod(x).second;}
    BigInteger BigInteger::pow(const long long& x) const {
        BigInteger res = 1, a = *this;
        for (long long t = x; t != 0; t >>= 1) {
            if (t & 1) res *= a;
            a *= a;
        }
        return res;
    }
    BigInteger BigInteger::pow(const long long& x, const BigInteger& p) const {
        BigInteger res = 1, a = *this % p;
        for (long long t = x; t != 0; t >>= 1) {
            if (t & 1) res = res * a % p;
            a = a * a % p;
        }
        return res;
    }

    BigInteger BigInteger::root(const long long& m) const {
        if (*this == 0 || m == 1) return *this;
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
        return x;
    }

    BigInteger BigInteger::gcd(const BigInteger& x) const {
        BigInteger a = *this, b = x;
        if (a < b) std::swap(a, b);
        if (b == 0) return a;
        int t = 0;
        while (a % 2 == 0 && b % 2 == 0) a = a.div2(), b = b.div2(), t++;
        while (b > 0) {
            if (a % 2 == 0) a = a.div2();
            else if (b % 2 == 0) b = b.div2();
            else a -= b;
            if (a < b) std::swap(a, b);
        }
        while (t--) a += a;
        return a;
    }
    BigInteger BigInteger::lcm(const BigInteger& x) const {return *this / gcd(x) * x;}

    BigInteger& BigInteger::operator+= (const BigInteger& x) {return *this = *this + x;}
    BigInteger& BigInteger::operator-= (const BigInteger& x) {return *this = *this - x;}
    BigInteger& BigInteger::operator*= (const BigInteger& x) {return *this = *this * x;}
    BigInteger& BigInteger::operator/= (const long long& x) {return *this = *this / x;}
    BigInteger& BigInteger::operator/= (const BigInteger& x) {return *this = *this / x;}
    BigInteger& BigInteger::operator%= (const long long& x) {return *this = *this / x;}
    BigInteger& BigInteger::operator%= (const BigInteger& x) {return *this = *this % x;}

    BigInteger BigInteger::operator<< (const long long& x) {
        if (x <= 0) return *this;
        BigInteger res = *this;
        for (long long i = 1; i <= x; i++) res += res;
        return res;
    }
    BigInteger BigInteger::operator>> (const long long& x) {
        if (x <= 0) return *this;
        BigInteger res = *this;
        for (long long i = 1; i <= x; i++) res = res.div2();
        return res;
    }
    BigInteger& BigInteger::operator<<= (const long long& x) {return *this = *this << x;}
    BigInteger& BigInteger::operator>>= (const long long& x) {return *this = *this >> x;}

    template <class F>
    inline BigInteger BigInteger::binary_op_helper(const BigInteger& x, const BigInteger& y, const F& func) {
        auto to_bin = [](BigInteger x) -> std::vector<bool> {
            if (x == 0) return {0};
            std::vector<bool> res;
            for (; x != 0; x = x.div2()) res.emplace_back(x.digits[1] & 1);
            return res;
        };
        std::vector<bool> a = to_bin(x), b = to_bin(y);
        int n = a.size(), m = b.size(), lim = std::max(n, m);
        std::vector<bool> res(lim, 0);
        for (int i = lim - 1; i >= 0; i--) 
            res[i] = func(i < n ? a[i] : 0, i < m ? b[i] : 0);
        std::reverse(res.begin(), res.end());
        return res;
    }
    BigInteger BigInteger::operator& (const BigInteger& x) {return binary_op_helper(*this, x, [](bool a, bool b) -> bool {return a & b;});}
    BigInteger BigInteger::operator| (const BigInteger& x) {return binary_op_helper(*this, x, [](bool a, bool b) -> bool {return a | b;});}
    BigInteger BigInteger::operator^ (const BigInteger& x) {return binary_op_helper(*this, x, [](bool a, bool b) -> bool {return a ^ b;});}
    BigInteger& BigInteger::operator&= (const BigInteger& x) {return *this = *this & x;}
    BigInteger& BigInteger::operator|= (const BigInteger& x) {return *this = *this | x;}
    BigInteger& BigInteger::operator^= (const BigInteger& x) {return *this = *this ^ x;}

    BigInteger& BigInteger::operator++ () {return *this += 1;}
    BigInteger BigInteger::operator++ (int) {BigInteger t = *this; return *this += 1, t;}
    BigInteger& BigInteger::operator-- () {return *this -= 1;}
    BigInteger BigInteger::operator-- (int) {BigInteger t = *this; return *this -= 1, t;}
} using namespace BGI;
```
### exgcd.hpp

```hpp
ll exgcd(ll a, ll b, ll &x, ll &y){
    if (b == 0){
        x = 1, y = 0;
        return a;
    }
    ll d = exgcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
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
    bool operator<(const coler_seg&a) const { return l < a.l; }
};
struct Chtholly : std::set<coler_seg> {
public:
    void add(int l, int r, DAT val) {
        auto itr = split(r + 1), itl = split(l);
        for (auto it = itl; it != itr; ++it) {
            it->val += val;
        }
    }
    void assign(int l, int r, DAT val){
        auto itr = split(r + 1), itl = split(l);
        erase(itl, itr);
        emplace(l, r, val);
    }
    ll kth(int l, int r, int rk) {
        auto itr = split(r + 1), itl = split(l);
        vector<pair<ll, int>> v;
        for (auto it = itl; it != itr; ++it) {
            v.emplace_back(it->val, it->r - it->l + 1);
        }
        MEION::sort(v);
        for (const auto &[val, sz] : v) {
            if (rk <= sz) return val;
            rk -= sz;
        }
        return LLMAX;
    }
    ll quis(int l, int r, int T, int mod) {
        auto itr = split(r + 1), itl = split(l);
        ll res = 0;
        for (auto it = itl; it != itr; ++it) {
            res = (res + (it->r - it->l + 1ll) * ksm((it->val) % mod, T, mod)) % mod;
        }
        return res;
    }
private:
    ll ksm(int a, int b, int mod) { ll res = 1; while (b) { if (b & 1) res = (res * a) % mod; a = 1ll * a * a % mod; b >>= 1; } return res % mod; }
    iterator split(int pos) {
        auto it = lower_bound(coler_seg(pos));
        if (it != end() and it->l == pos) return it;
        coler_seg tmp = *--it;
        erase(it);
        emplace(tmp.l, pos - 1, tmp.val);
        return emplace(pos, tmp.r, tmp.val).first;
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
        return cnt;
    }
    int dir(int x) { return ch[fa[x]][1] == x; }
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
            return;
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
            else if (x <= siz[ch[o][0]] + num[o]) return val[o];
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
    int operator[](int x) { return ff(x); }
    int size(int x) { return sz[ff(x)]; }
    bool merge(int x, int y) { 
        x = ff(x), y = ff(y); 
        if (x == y) return false; 
        --comp; 
        sz[x] += sz[y], sz[y] = 0; fa[y] = x; 
        return true; 
    }
private:
    int n, comp;
    std::vector<int> fa, sz;
    int ff(int x) { 
        while (x != fa[x]) x = fa[x] = fa[fa[x]]; 
        return x; 
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
        return v[i >> 6] >> (i & 63) & 1;
    }
	bitarray operator &=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] &= b.v[i];
        }
        return *this;
	}
	bitarray operator |=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] |= b.v[i];
        }
        return *this;
	}
	bitarray operator ^=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] ^= b.v[i];
        }
        return *this;
	}
	bitarray operator &(const bitarray &b) {
        return bitarray(*this) &= b;
	}
	bitarray operator |(const bitarray &b) {
        return bitarray(*this) |= b;
	}
	bitarray operator ^(const bitarray &b) {
        return bitarray(*this) ^= b;
	}
    bitarray operator ~() const {
		bitarray ret(*this);
        for (int i = 0, ed = sz; i < ed; ++i) {
            ret.v[i] = ~ret.v[i];
        }
        return ret;
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
        return bitarray(*this) <<= t;
    }
    bitarray operator >>(const int t) {
        return bitarray(*this) >>= t;
    }
    std::string to_string() {
        std::string ans;
        for (int i = 0; i < N; ++i) {
            ans += '0' + (*this)[i];
        }
        return ans;
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
        return v[i >> 6] >> (i & 63) & 1;
    }
	bitvector operator &=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] &= b.v[i];
        }
        return *this;
	}
	bitvector operator |=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] |= b.v[i];
        }
        return *this;
	}
	bitvector operator ^=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] ^= b.v[i];
        }
        return *this;
	}
	bitvector operator &(const bitvector &b) {
        return bitvector(*this) &= b;
	}
	bitvector operator |(const bitvector &b) {
        return bitvector(*this) |= b;
	}
	bitvector operator ^(const bitvector &b) {
        return bitvector(*this) ^= b;
	}
    bitvector operator ~() const {
		bitvector ret(*this);
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            ret.v[i] = ~ret.v[i];
        }
        return ret;
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
        return bitvector(*this) <<= t;
    }
    bitvector operator >>(const int t) {
        return bitvector(*this) >>= t;
    }
    std::string to_string() {
        std::string ans;
        for (int i = 0; i < n; ++i) {
            ans += '0' + (*this)[i];
        }
        return ans;
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
        return history.size(); 
    }
    void rollback(int t) { 
        for (int i = time() - 1; i >= t; --i) { 
            auto& [idx, v] = history[i]; dat[idx] = v; 
        } 
        history.resize(t); 
    }
    T get(int idx) { 
        return dat[idx]; 
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
        return res; }
    T operator[](int idx) { 
        return dat[idx]; 
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
    int size() { return used.size() / 2 - cap; }
    int index(const ull &k) {
        int i = 0;
        for (i = hash(k); used[i] and key[i] != k; i = (i + 1) & msk) {}
        return i;
    }

    Val& operator[](const ull &k) {
        if (cap == 0) extend();
        int i = index(k);
        if (not used[i]) { used[i] = 1, key[i] = k, val[i] = Val{}, --cap; }
        return val[i];
    }

    Val get(const ull &k, Val default_value) {
        int i = index(k);
        return (used[i] ? val[i] : default_value);
    }

    bool count(const ull &k) {
        int i = index(k);
        return used[i] and key[i] == k;
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
        return (x ^ (x >> 31)) & msk;
    }

    void extend() {
        vector<pair<ull, Val>> dat;
        dat.reserve(used.size() / 2 - cap);
        for (int i = 0, ed = used.size(); i < ed; ++i) {
            if (used[i]) dat.emplace_back(key[i], val[i]);
        }
        build(dat.size() << 1);
        for (auto &[a, b] : dat) (*this)[a] = b;
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
        return v;
    }
    ll quis(ll v = 0ll) {
        for (int i = B; i--; ) {
            if (not (v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        return v;
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
        build(v.size(), [&](int i) -> T { return v[i]; });
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

    T sum_all() { return total; }
    T sum(int k) { return pre_sum(k); }
    T pre_sum(int k) {
        chmin(k, n);
        T res(0);
        for (; k > 0; k -= k & -k) {
            res += dat[k - 1];
        }
        return res;
    }
    T sum(int l, int r) {
        chmax(l, 0);
        chmin(r, n);
        if (l == 0) return pre_sum(r);
        T pos = T(0), neg = T(0);
        while (l < r) {
            pos += dat[r - 1];
            r -= r & -r;
        }
        while (r < l) {
            neg += dat[l - 1];
            l -= l & -l;
        }
        return pos - neg;
    }
    vector<T> get_all() {
        vector<T> res(n);
        for (int i = 0; i < n; ++i) {
            res[i] = sum(i, i + 1);
        }
        return res;
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

    int sum_all() { return bit.sum_all(); }
    int pre_sum(int k) {
        int ans = bit.sum(k / 64);
        ans += popcount(dat[k / 64] & ((1ull << (k % 64)) - 1));
        return ans;
    }
    int sum(int k) { return pre_sum(k); }
    int sum(int l, int r) {
        if (l == 0) return pre_sum(r);
        int ans = 0;
        ans -= popcount(dat[l / 64] & ((1ull << (l % 64)) - 1));
        ans += popcount(dat[r / 64] & ((1ull << (r % 64)) - 1));
        ans += bit.sum(l / 64, r / 64);
        return ans;
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
        return v;
    }
    int size(int v) { 
        return -dat.get((*this)[v]); 
    }
    int time() { 
        return dat.time(); 
    }
    void rollback(int t) { 
        dat.rollback(t); 
    }
    bool merge(int a, int b) {
        a = (*this)[a], b = (*this)[b];
        if (a == b) {
            return false;
        }
        if (dat.get(a) > dat.get(b)) {
            std::swap(a, b);
        }
        dat.set(a, dat.get(a) + dat.get(b));
        dat.set(b, a);
        return true;
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
	T top() { return p.top(); }
	bool empty() { return p.empty(); }
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
        auto [a, b] = dat[k >> 5];
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
            for (auto&& i : I) {
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
    int count(int L, int R, T a, T b, T xor_val = 0) { return prefix_count(L, R, b, xor_val) - prefix_count(L, R, a, xor_val); }
    // xor した結果で [0, x) に収まるものを数える
    int prefix_count(int L, int R, T x, T xor_val = 0) {
        if (xor_val != 0) assert(set_log);
        x = (COMPRESS ? std::distance((key).begin(), std::lower_bound(key.begin(), key.end(), (x))) : x);
        if (x >= (1 << lg)) return R - L;
        int ret = 0;
        for (int d = lg - 1; d >= 0; --d) {
            bool add = (x >> d) & 1;
            bool f = ((x ^ xor_val) >> d) & 1;
            if (add) ret += bv[d].rank(R, !f) - bv[d].rank(L, !f);
            L = bv[d].rank(L, f) + (f ? mid[d] : 0);
            R = bv[d].rank(R, f) + (f ? mid[d] : 0);
        }
        return ret;
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
        return (COMPRESS ? key[ret] : ret);
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
            if (id[i << 1] == id[i << 1 | 1]) return false;
            if (id[i << 1] < id[i << 1 | 1]) {
                ans[i] = 1;
            }
        }
        return true;
    }
    bitvector answer() { return ans; }
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
            auto cmp = [&] (int i, int k) {
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
        if (n < 1) return false;
        for (int i = 2, ed = n; i * i <= ed; ++i) {
            if (n % i == 0) return false;
        }
        return true;
    }
    int guidingstar_find_pr(int n) {
        while (not guidingstar_ckpr(n)) ++n;
        return n;
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
        return { (h[r].first + 1ll * (getmod::m1 - h[l].first) * p[r - l].first) % getmod::m1,
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
                if (at(p).len + 1 == at(q).len) return q;
                int cp = size();
                push_back(at(q));
                back().len = at(p).len + 1;
                while (~p and at(p)[c] == q) {
                    at(p)[c] = cp;
                    p = at(p).link;
                }
                at(q).link = cp;
                return cp;
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
            return pla;
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
            return pla;
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
        return {x + p.x, y + p.y};
    }
    point operator-(point p) const {
        return {x - p.x, y - p.y};
    }
    bool operator==(point p) const {
        return x == p.x and y == p.y;
    }
    bool operator!=(point p) const {
        return x != p.x or y != p.y;
    }
    point operator-() const {
        return {-x, -y};
    }
    point operator*(T t) const {
        return {x * t, y * t};
    }
    point operator/(T t) const {
        return {x / t, y / t};
    }
 
    bool operator<(point p) const {
        if (x != p.x) return x < p.x;
        return y < p.y;
    }
    bool operator>(point p) const {
        if (x != p.x) return x > p.x;
        return y > p.y;
    }
    T dot(const point &other) const {
        return x * other.x + y * other.y;
    }
    T det(const point &other) const {
        return x * other.y - y * other.x;
    }
    T square() const {
        return x * x + y * y;
    }
 
    RE length() { return sqrtl(x * x + y * y); }
    RE angle() { return std::atan2(y, x); }
 
    point rotate(double theta) {
        static_assert(not std::is_integral<T>::value);
        RE c = std::cos(theta), s = std::sin(theta);
        return point{c * x - s * y, s * x + c * y};
    }
    point rot90(bool ccw = 1) {
        return (ccw ? point{-y, x} : point{y, -x});
    }
};
 
template <typename T>
std::istream& operator>>(std::istream& is, point<T>& any) {
    is >> any.x >> any.y;
    return is;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const point<T>& any) {
    os << any.x << ' ' << any.y;
    return os;
}

// A -> B -> Cと進むときに，左转为 +1，右转为 -1。
template<typename T>
int ccw(point<T> a, point<T> b, point<T> c) {
    T x = (b - a).det(c - a);
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

template <typename REAL = long double, typename T, typename U>
REAL dist(point<T> a, point<U> b) {
    REAL dx = REAL(a.x) - REAL(b.x);
    REAL dy = REAL(a.y) - REAL(b.y);
    return std::sqrt(dx * dx + dy * dy);
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
        return a * p.y + b * p.y + c;
    }

    template <typename U>
    T eval(U x, U y) const {
        return a + x + b * y + c;
    }

    void normalize() {
        static_assert(std::is_same_v<T, int> or std::is_same_v<T, long long>);
        T gcd = std::gcd(std::gcd(std::abs(a), std::abs(b)), std::abs(c));
        a /= gcd, b /= gcd, c /= gcd;
    }

    bool parallel(line other) const {
        return a * other.b - b * other.a == 0; 
    }
    bool is_orthoginal(line other) const {
        return a * other.a + b * other.b == 0;
    }
};

template <typename T>
struct segment {
    point<T> a, b;
    
    segment(point<T> a, point<T> b) : a(a), b(b) {}
    segment(T x1, T y1, T x2, T y2) : segment(point<T>(x1, y1), point<T>(x2, y2)) {}

    bool contain(point<T> c) const {
        T det = (c - a).det(b - a);
        if (det != 0) return 0;
        return (c - a).dot(b - a) >= 0 and (c - b).dot(a - b) >= 0;
    }

    line<T> to_line() { return line(a, b); }
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
        return dx * dx + dy * dy <= r * r;
    }
};

// 不平行仮定
template <typename REAL = long double, typename T>
point<REAL> cross_point(const line<T> l1, const line<T> l2) {
    T det = l1.a * l2.b - l1.b * l2.a;
    assert(det != 0);
    REAL x = -REAL(l1.c) * l2.b + REAL(l1.b) * l2.c;
    REAL y = -REAL(l1.a) * l2.c + REAL(l1.c) * l2.a;
    return point<REAL>(x / det, y / det);
}
template <typename REAL = long double, typename T>
point<REAL> line_x_line(const line<T> l1, const line<T> l2) {
    return cross_point<REAL, T>(l1, l2);
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
        if (l1.eval(s2.a) != 0) return 0;
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
        if (a < b) return 2;
        if (a > b) return 0;
        return (include_ends ? 1 : 0);
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
    return (ok1 and ok2 ? 1 : 0);
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
    if (D < 0) return {};
    REAL sqD = sqrtl(D);
    REAL y1 = (-2 * b * c + sqD) / (2 * (a * a + b * b));
    REAL y2 = (-2 * b * c - sqD) / (2 * (a * a + b * b));
    REAL x1 = (-b * y1 - c) / a;
    REAL x2 = (-b * y2 - c) / a;
    if (sw) std::swap(x1, y1), std::swap(x2, y2);
    x1 += C.O.x, x2 += C.O.x;
    y1 += C.O.y, y2 += C.O.y;
    if (D == 0) {
        return {point<REAL>(x1, y1)};
    }
    return {point<REAL>(x1, y1), point<REAL>(x2, y2)};
}

template <typename REAL, typename T>
std::tuple<bool, point<T>, point<T>> cross_point_circle(circle<T> C1,
                                                        circle<T> C2) {
    using P = point<T>;
    P O {0, 0};
    P A = C1.O, B = C2.O;
    if (A == B) return {false, O, O};
    T d = (B - A).length();
    REAL cos_val = (C1.r * C1.r + d * d - C2.r * C2.r) / (2 * C1.r * d);
    if (cos_val < -1 || 1 < cos_val) return {false, O, O};
    REAL t = std::acos(cos_val);
    REAL u = (B - A).angle();
    P X = A + P {C1.r * std::cos(u + t), C1.r * std::sin(u + t)};
    P Y = A + P {C1.r * std::cos(u - t), C1.r * std::sin(u - t)};
    return {true, X, Y};
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
    return circle<REAL>(O.x, O.y, r);
}
```
### 3-angle_sort.hpp

```hpp
#pragma once
#include "1-base.hpp"

template <typename T>
int lower_or_upper(const point<T> &p) {
    if (p.y != 0) return (p.y > 0 ? 1 : -1);
    if (p.x > 0) return -1;
    if (p.x < 0) return 1;
    return 0;
}

template <typename T>
int angle_cmp_3(const point<T> &L, const point<T> &R) {
    int a = lower_or_upper(L), b = lower_or_upper(R);
    if (a != b) return (a < b ? -1 : +1);
    T det = L.det(R);
    if (det > 0) return -1;
    if (det < 0) return 1;
    return 0;
}

template <typename T>
vector<int> angle_sort(const vector<point<T>> &v) {
    vector<int> rk(v.size());
    std::iota(rk.begin(), rk.end(), 0);
    sort(rk, [&](auto &L, auto &R) -> bool{
        return (angle_cmp_3(v[L], v[R]) == -1);
    });
    return rk;
}

template <typename T>
vector<int> angle_sort(const vector<pair<T, T>> &v) {
    vector<point<T>> tmp(v.size());
    for (int i = 0, ed = v.size(); i < ed; ++i) {
        tmp[i] = point<T>(v[i]);
    }
    return angle_sort(tmp);
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

    auto square = [&] (int i, int j) -> T {
        return (points[j] - points[i]).square();
    };

    T best = square(0, 1);
    pair<int, int> res(0, 1);
    T w =  sqrtl(best);
    
    vector<int> nxt(n, -1);
    
    auto ins = [&] (int i) -> void {
        ull k = hash_pair(pair{points[i].x / w, points[i].y / w});
        nxt[i] = ma.get(k, -1);
        ma[k] = i;
    };

    auto quis = [&] (int i) -> bool {
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
        return upd;
    };

    if (best == T(0)) {
        res.first = id[res.first], res.second = id[res.second];
        return res;
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
    return res;
}

template <typename T>
pair<int, int> closest_pair2(vector<point<T>> points) {
    using RE = long double;
    const int n = points.size();
    if (n == 1) return pair(0, 0);
    ld rd = MeIoN_random_hash::rng(114514) % 360 * 0.114514;
    ld SIN = std::cos(rd), COS = std::sin(rd);
    vector<int> id(n);
    for (int i = 0; i < n; ++i) id[i] = i;
    sort(id, [&](auto &a, auto &b) -> bool {
        return points[a].x * COS - points[a].y * SIN < 
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
    return ans;
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
    if (n == 1) return {0};
    if (n == 2) {
        if (p[0] < p[1]) return {0, 1};
        if (p[0] > p[1]) return {1, 0};
        return {0};
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
    auto check = [&](int i, int j, int k) -> bool {
        T det = (p[j] - p[i]).det(p[k] - p[i]);
        if constexpr (allow_180) {
            return det >= 0;
        }
        return det > T(0);
    };
    auto cal = [&]() {
        vector<int> res;
        for (const auto &k : id) {
            while (res.size() > 1) {
                auto i = res.end()[-2];
                auto j = res.end()[-1];
                if (check(i, j, k)) break;
                res.pop_back();
            }
            res.emplace_back(k);
        }
        return res;
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
    return res;
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
        return m % n;
    }

    int nxt_idx(int i) { return (i + 1 == n ? 0 : i + 1); }
    int pre_idx(int i) { return (i == 0 ? n - 1 : i - 1); }

    // 中：1, 境界：0, 外：-1. test した.
    int side(P p) {
        int l = 1, r = n - 1;
        T a = (points[l] - points[0]).det(p - points[0]);
        T b = (points[r] - points[0]).det(p - points[0]);
        if (a < 0 or b > 0) return -1;
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
        if (x < 0) return -1;
        if (x > 0) return 1;
        // on triangle p[0] p[L] p[R]
        if (p == points[0]) return 0;
        if (c != 0 and a == 0 and l != 1) return 1;
        if (c != 0 and b == 0 and r != n - 1) return 1;
        return 0;
    }

    // return {min, idx} 点积最小值 O(log)
    pair<T, int> min_dot(P p) {
        int idx = periodic_min_comp([&](int i, int k) -> bool {
            return points[i].dot(p) < points[k].dot(p);
        });
        return {points[idx].dot(p), idx};
    }

    // return {max, idx} 点积最大值 O(log)
    pair<T, int> max_dot(P p) {
        int idx = periodic_min_comp([&](int i, int k) -> bool {
            return points[i].dot(p) > points[k].dot(p);
        });
        return {points[idx].dot(p), idx};
    }

    // 计算从一个点 p 观看多边形时，可以看到的多边形的范围
    // 该函数返回两个索引，表示可以看到的范围（考虑反向偏角）
    pair<int, int> visible_range(P p) {
        int a = periodic_min_comp([&](int i, int k) -> bool {
            return ((points[i] - p).det(points[k] - p) < 0);
        });
        int b = periodic_min_comp([&](int i, int k) -> bool {
            return ((points[i] - p).det(points[k] - p) > 0);
        });
        if ((p - points[a]).det(p - points[pre_idx(a)]) == T(0)) a = pre_idx(a);
        if ((p - points[b]).det(p - points[nxt_idx(b)]) == T(0)) b = nxt_idx(b);
        return {a, b};
    }

    // 线段是否与凸多边形相交
    bool check_cross(P pa, P pb) {
        for (int i = 0; i < 2; ++i) {
            std::swap(pa, pb);
            auto [a, b] = visible_range(pa);
            if ((points[a] - pa).det(pb - pa) >= 0) return false;
            if ((points[b] - pa).det(pb - pa) <= 0) return false;
        }
        return true;
    }

    vector<T> AREA;

    // point[i,...,j] (inclusive) 面积
    T area_between(int i, int k) {
        assert(-1 < i and i < n);
        assert(-1 < k and k < n);
        if (i > k) k += n;
        if (AREA.empty()) build_AREA();
        return AREA[k] - AREA[i] + (points[k % n].det(points[i]));
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
        for (const auto p : A)
            assert(std::max(std::abs(p.x), std::abs(p.y)) < limit);
        for (const auto p : B)
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
        if (d == 0) return 0;
        if (d > 0) {
            return tri[i][j] + tri[j][k] - tri[i][k] - seg[i][k];
        }
        int x = tri[i][k] - tri[i][j] - tri[j][k];
        return x - seg[i][j] - seg[j][k] - points[j];
    }

    int count2(int i, int j) {
        i = new_idx[i], j = new_idx[j];
        if (i > j) std::swap(i, j);
        return seg[i][j];
    }

   private:
    P take_origin() {
        // 不要让OAiAj和OAiBj在同一直线上
        // fail prob: at most N(N+M)/LIM
        // return P {-limit, MeIoN_random_hash::rng_64(-limit, limit)};
        return P {-limit, MeIoN_random_hash::rng(-limit, limit)};
    }

    void build() {
        P O = take_origin();
        for (auto &p : A) {
            p = p - O;
        }
        for (auto &p : B) {
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
                auto x = (ng + ok) >> 1;
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
                 [&](auto &a, auto &b) -> bool { return C[a].det(C[b]) > 0; });
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
                        return(n == 0 ? true : C[n - 1].det(p) > 0);
                    }, 0, m + 1);
                int ub = binary_search(
                    [&](int n) -> bool {
                        return(n == 0 ? true : C[n - 1].det(p) >= 0);
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
    return std::sqrt(dx * dx + dy * dy);
}

template <typename REAL = ld, typename T, typename U>
REAL distance(segment<T> S, point<U> P) {
    point<T> A = S.a, B = S.b;
    bool b1 = (B - A).dot(P - A) >= 0;
    bool b2 = (A - B).dot(P - B) >= 0;
    if (b1 and not b2) {
        return distance<REAL, T, T>(B, P);
    }
    if (not b1 and b2) {
        return distance<REAL, T, T>(A, P);
    }
    line<T> L = S.to_line();
    return REAL(std::abs(L.eval(P))) / std::sqrt(REAL(L.a) * L.a + REAL(L.b) * L.b);
}

template <typename REAL, typename T>
REAL distance(segment<T> s1, segment<T> s2) {
    if (count_cross<T>(s1, s2, true)) return REAL(0);
    REAL res = distance<REAL, T, T>(s1, s2.a);
    chmin(res, distance<REAL, T, T>(s1, s2.b));
    chmin(res, distance<REAL, T, T>(s2, s1.a));
    chmin(res, distance<REAL, T, T>(s2, s1.b));
    return res;
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

    auto upd = [&](int i, int j) -> void {
        point<T> p = points[i] - points[j];
        ll d = p.dot(p);
        if (chmax(best, d)) ans = pair(i, j);
    };
    upd(0, 1);

    vector<int> id = convex_hull(points);
    int n = id.size();
    if (n == 1) return ans;
    if (n == 2) return pair(id[0], id[1]);
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
    return ans;
}
```
### 10-triangle_area.hpp

```hpp
#pragma once
#include "1-base.hpp"

template <typename REAL = long double, typename T>
REAL triangle_area(point<T> a, point<T> b, point<T> c) {
    return std::abs((b - a).det(c - a)) * 0.5L;
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
    return Circle<REAL>(x, y, r);
}
```

<div style="page-break-after: always;"></div>


## random

### random.hpp

```hpp
namespace MeIoN_random_hash {
    std::mt19937 RNG(std::chrono::steady_clock::now().time_since_epoch().count());
    uint rng(uint limit) { return RNG() % limit; }
    int rng(int l, int r) { return l + RNG() % (r - l); }
    std::mt19937_64 RNG_64(std::chrono::steady_clock::now().time_since_epoch().count());
    ull rng_64(ull limit) { return RNG_64() % limit; }
    ll rng_64(ll l, ll r) { return l + RNG_64() % (r - l); }

    using m1 = modint<998244353>;
    using m2 = modint<1000000007>;

    namespace get_prim {

        constexpr ull md = (1ull << 61) - 1;

        static inline constexpr ull modmul(const ull &a, const ull &b) {
            u128 d = u128(a) * b;
            ull ret = (ull(d) & md) + ull(d >> 61);
            return ret >= md ? ret - md : ret;
        }
        
        static ull modpow(ull a, ull b) {
            ull r = 1;
            for (a %= md; b; a = modmul(a, a), b >>= 1) r = modmul(r, a);
            return r;
        }

        static bool is_primitive(ull x) {
            for (auto &d :
                vector<ull> {2, 3, 5, 7, 11, 13, 31, 41, 61, 151, 331, 1321})
                if (modpow(x, (md - 1) / d) <= 1) return false;
            return true;
        }

        static ull get_basis() {
            static auto rand_time =
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::high_resolution_clock::now().time_since_epoch())
                    .count();
            static std::mt19937_64 rng(rand_time);
            ull ret;
            while (is_primitive(ret = rng() % (md - 1) + 1) == false);
            return ret;
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
        for (auto &[x, y] : v) {
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
                auto [a, b] = cand[i];
                v.emplace_back(a, b);
                break;
            }
        }
        random_relabel(n, v);
        return v;
    }

    template <typename T>
    ull hash_pair(const pair<T, T> &X) {
        static ll hash_base = RNG_64();
        if (hash_base == 0) hash_base = RNG_64();
        return hash_base * X.first + X.second;
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
        return pair(h1.val, h2.val);
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
        return pair(h1.val, h2.val);
    }

    // https://uoj.ac/problem/763
    struct rooted_tree_hash {
        vector<vector<int>> v;
        int n;
        vector<ull> hash;
        vector<int> dis;

        static vector<ull> &xs() {
            static vector<ull> _xs;
            return _xs;
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
            return dis[n] = dp;
        }
    };
}
```

<div style="page-break-after: always;"></div>


## mod

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
    if (int(c.val) < n) return f[c.val];
    auto a = f;
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
    return res;
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
        return res;
    }
};
```
### modint.hpp

```hpp
template <int mod>
struct modint {
    static constexpr bool is_mod_int = true;
    static constexpr unsigned umod = unsigned(mod);
    static_assert(umod < unsigned(1) << 31);
    int val;
    static modint raw(unsigned v) {
        modint x;
        x.val = v;
        return x;
    }
    constexpr modint(const ll val = 0) noexcept : val(val >= 0 ? val % mod : (mod - (-val) % mod) % mod) { }
    bool operator<(const modint& other) const { return val < other.val; }
    modint& operator+=(const modint& p) {
        if ((val += p.val) >= mod)
            val -= mod;
        return* this;
    }
    modint& operator-=(const modint& p) {
        if ((val += mod - p.val) >= mod)
            val -= mod;
        return* this;
    }
    modint& operator*=(const modint& p) {
        val = (int)(1LL * val * p.val % mod);
        return* this;
    }
    modint& operator/=(const modint& p) {
        *this *= p.inv();
        return* this;
    }
    modint operator-() const { return modint::raw(val ? mod - val : unsigned(0)); }
    modint operator+(const modint& p) const { return modint(*this) += p; }
    modint operator-(const modint& p) const { return modint(*this) -= p; }
    modint operator*(const modint& p) const { return modint(*this) *= p; }
    modint operator/(const modint& p) const { return modint(*this) /= p; }
    bool operator==(const modint& p) const { return val == p.val; }
    bool operator!=(const modint& p) const { return val != p.val; }
    friend std::istream& operator>>(std::istream& is, modint& p) {
        ll x;
        is >> x;
        p = x;
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, modint p) { return os << p.val; }
    modint inv() const {
        int a = val, b = mod, u = 1, v = 0, t;
        while (b > 0)
            t = a / b, std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        return modint(u);
    }
    modint ksm(ll n) const {
        modint ret(1), mul(val);
        while (n > 0) {
            if (n & 1)
                ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }
    static constexpr int get_mod() { return mod; }
    static constexpr pair<int, int> ntt_info() {
        if (mod == 120586241) return {20, 74066978};
        if (mod == 167772161) return {25, 17};
        if (mod == 469762049) return {26, 30};
        if (mod == 754974721) return {24, 362};
        if (mod == 880803841) return {23, 211};
        if (mod == 943718401) return {22, 663003469};
        if (mod == 998244353) return {23, 31};
        if (mod == 1004535809) return {21, 836905998};
        if (mod == 1045430273) return {20, 363};
        if (mod == 1051721729) return {20, 330};
        if (mod == 1053818881) return {20, 2789};
        return { -1, -1 };
    }
    static constexpr bool can_ntt() { return ntt_info().first != -1; }
};
ll mod_inv(ll val, ll mod) {
    if (mod == 0) return 0;
    mod = std::abs(mod);
    val %= mod;
    if (val < 0) val += mod;
    ll a = val, b = mod, u = 1, v = 0, t;
    while (b > 0) {
        t = a / b;
        std::swap(a -= t * b, b), std::swap(u -= t * v, v);
    }
    if (u < 0) u += mod;
    return u;
}
constexpr unsigned mod_pow_constexpr(ull a, ull n, unsigned mod) {
    a %= mod;
    ull res = 1;
    for (int _ = 0; _ < 32; ++_) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod, n /= 2;
    }
    return res;
}
template <typename T, unsigned p0, unsigned p1, unsigned p2>
T CRT3(ull a0, ull a1, ull a2) {
    static_assert(p0 < p1 && p1 < p2);
    static constexpr ull x0_1 = mod_pow_constexpr(p0, p1 - 2, p1);
    static constexpr ull x01_2 = mod_pow_constexpr(ull(p0) * p1 % p2, p2 - 2, p2);
    ull c = (a1 - a0 + p1) * x0_1 % p1;
    ull a = a0 + c * p0;
    c = (a2 - a % p2 + p2) * x01_2 % p2;
    return T(a) + T(c) * T(p0) * T(p1);
}
```

<div style="page-break-after: always;"></div>
