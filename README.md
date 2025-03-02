#### Categorization
![alt text](Z_some_tools/cover.jpg)
### MeIoN's XCPC_Library
#### qq : 60422310
#### Codeforces : MeIoN_is_UMP45, white_unicorn
#### Luogu : MeIoN
#### Atcoder : MeIoN
##### using cpp_merge to bind cpp file (https://github.com/FastAlien/cpp-merge)


# 目录

- [Z\_H/MeIoN\_H.hpp](#z_hmeion_hhpp)
- [Z\_H/MeIoN\_IO.hpp](#z_hmeion_iohpp)
- [Z\_H/MeIoN\_PRET.hpp](#z_hmeion_prethpp)
- [Z\_H/MeIoN\_debug.hpp](#z_hmeion_debughpp)
- [Z\_H/fast\_io.hpp](#z_hfast_iohpp)
- [ds/LinearBasis.hpp](#dslinearbasishpp)
- [ds/Wavelet\_Matrix.hpp](#dswavelet_matrixhpp)
- [ds/a\_monoid/max\_add.hpp](#dsa_monoidmax_addhpp)
- [ds/a\_monoid/max\_assign.hpp](#dsa_monoidmax_assignhpp)
- [ds/a\_monoid/min\_add.hpp](#dsa_monoidmin_addhpp)
- [ds/a\_monoid/minidx\_add.hpp](#dsa_monoidminidx_addhpp)
- [ds/a\_monoid/minmax\_add.hpp](#dsa_monoidminmax_addhpp)
- [ds/a\_monoid/minmax\_rev.hpp](#dsa_monoidminmax_revhpp)
- [ds/a\_monoid/minmincnt\_add.hpp](#dsa_monoidminmincnt_addhpp)
- [ds/a\_monoid/sum\_add.hpp](#dsa_monoidsum_addhpp)
- [ds/a\_monoid/sum\_assign.hpp](#dsa_monoidsum_assignhpp)
- [ds/a\_monoid/sum\_mul.hpp](#dsa_monoidsum_mulhpp)
- [ds/a\_monoid/sum\_rev.hpp](#dsa_monoidsum_revhpp)
- [ds/bit\_vec.hpp](#dsbit_vechpp)
- [ds/chtholly.hpp](#dschthollyhpp)
- [ds/dsu.hpp](#dsdsuhpp)
- [ds/fenw.hpp](#dsfenwhpp)
- [ds/hashmap.hpp](#dshashmaphpp)
- [ds/heap.hpp](#dsheaphpp)
- [ds/monoid/add.hpp](#dsmonoidaddhpp)
- [ds/monoid/add\_array.hpp](#dsmonoidadd_arrayhpp)
- [ds/monoid/add\_pair.hpp](#dsmonoidadd_pairhpp)
- [ds/monoid/and.hpp](#dsmonoidandhpp)
- [ds/monoid/assign.hpp](#dsmonoidassignhpp)
- [ds/monoid/gcd.hpp](#dsmonoidgcdhpp)
- [ds/monoid/max.hpp](#dsmonoidmaxhpp)
- [ds/monoid/max\_idx.hpp](#dsmonoidmax_idxhpp)
- [ds/monoid/min.hpp](#dsmonoidminhpp)
- [ds/monoid/min\_idx.hpp](#dsmonoidmin_idxhpp)
- [ds/monoid/minmax.hpp](#dsmonoidminmaxhpp)
- [ds/monoid/minmincnt.hpp](#dsmonoidminmincnthpp)
- [ds/monoid/mul.hpp](#dsmonoidmulhpp)
- [ds/monoid/reverse.hpp](#dsmonoidreversehpp)
- [ds/monoid/sum.hpp](#dsmonoidsumhpp)
- [ds/monoid/xor.hpp](#dsmonoidxorhpp)
- [ds/range\_mex\_query.hpp](#dsrange_mex_queryhpp)
- [ds/rectangle\_union.hpp](#dsrectangle_unionhpp)
- [ds/rollback\_array.hpp](#dsrollback_arrayhpp)
- [ds/rollback\_dsu.hpp](#dsrollback_dsuhpp)
- [ds/seg/dynamic\_seg.hpp](#dssegdynamic_seghpp)
- [ds/seg/lazy\_dynamic\_seg.hpp](#dsseglazy_dynamic_seghpp)
- [ds/seg/lazy\_seg\_base.hpp](#dsseglazy_seg_basehpp)
- [ds/seg/seg\_base.hpp](#dssegseg_basehpp)
- [ds/sparse\_table/st.hpp](#dssparse_tablesthpp)
- [ds/splay.hpp](#dssplayhpp)
- [ds/sqrt\_tree.hpp](#dssqrt_treehpp)
- [ds/st\_table.hpp](#dsst_tablehpp)
- [flow/hungarian.hpp](#flowhungarianhpp)
- [flow/max\_flow.hpp](#flowmax_flowhpp)
- [flow/max\_flow\_min\_cost.hpp](#flowmax_flow_min_costhpp)
- [flow/maxflow.hpp](#flowmaxflowhpp)
- [geo/1-base.hpp](#geo1-basehpp)
- [geo/10-triangle\_area.hpp](#geo10-triangle_areahpp)
- [geo/11-in\_circle.hpp](#geo11-in_circlehpp)
- [geo/12-line\_inside\_polygon.hpp](#geo12-line_inside_polygonhpp)
- [geo/13-manhattan\_mst.hpp](#geo13-manhattan_msthpp)
- [geo/14-max\_norm\_sum.hpp](#geo14-max_norm_sumhpp)
- [geo/15-minkowski\_sum.hpp](#geo15-minkowski_sumhpp)
- [geo/16-out\_circle.hpp](#geo16-out_circlehpp)
- [geo/17-minimum\_enclosing\_circle.hpp](#geo17-minimum_enclosing_circlehpp)
- [geo/2-apollonian\_circle.hpp](#geo2-apollonian_circlehpp)
- [geo/3-angle\_sort.hpp](#geo3-angle_sorthpp)
- [geo/4-closest\_pair.hpp](#geo4-closest_pairhpp)
- [geo/5-hull.hpp](#geo5-hullhpp)
- [geo/6-convex\_polygon.hpp](#geo6-convex_polygonhpp)
- [geo/7-points\_in\_triangles.hpp](#geo7-points_in_triangleshpp)
- [geo/8-distance.hpp](#geo8-distancehpp)
- [geo/9-furthest\_pair.hpp](#geo9-furthest_pairhpp)
- [graph/Apck/Basic.hpp](#graphapckbasichpp)
- [graph/Apck/bellman\_ford.hpp](#graphapckbellman_fordhpp)
- [graph/Apck/block\_cut\_tree.hpp](#graphapckblock_cut_treehpp)
- [graph/Apck/dag\_path\_cover.hpp](#graphapckdag_path_coverhpp)
- [graph/Apck/dijkstra.hpp](#graphapckdijkstrahpp)
- [graph/Apck/dominator\_tree.hpp](#graphapckdominator_treehpp)
- [graph/Apck/find\_centroids.hpp](#graphapckfind_centroidshpp)
- [graph/Tree/Basic.hpp](#graphtreebasichpp)
- [graph/Tree/dsu\_on\_tree.hpp](#graphtreedsu_on_treehpp)
- [graph/Tree/fast\_lca.hpp](#graphtreefast_lcahpp)
- [graph/Tree/tree\_monoid.hpp](#graphtreetree_monoidhpp)
- [graph/Tree/tree\_monoid\_lazy.hpp](#graphtreetree_monoid_lazyhpp)
- [graph/bellman\_ford.hpp](#graphbellman_fordhpp)
- [graph/dijkstra.hpp](#graphdijkstrahpp)
- [graph/find\_cycle\_directed.hpp](#graphfind_cycle_directedhpp)
- [graph/floyd.hpp](#graphfloydhpp)
- [graph/scc.hpp](#graphscchpp)
- [graph/triangle\_counting.hpp](#graphtriangle_countinghpp)
- [graph/two\_sat.hpp](#graphtwo_sathpp)
- [math/Big\_int.hpp](#mathbig_inthpp)
- [math/counting/count\_rectangle.hpp](#mathcountingcount_rectanglehpp)
- [math/crt.hpp](#mathcrthpp)
- [math/exgcd.hpp](#mathexgcdhpp)
- [math/line/transpose.hpp](#mathlinetransposehpp)
- [math/line/vector\_space.hpp](#mathlinevector_spacehpp)
- [math/mat.hpp](#mathmathpp)
- [math/mod/barrett.hpp](#mathmodbarretthpp)
- [math/mod/count\_terms.hpp](#mathmodcount_termshpp)
- [math/mod/differentiate.hpp](#mathmoddifferentiatehpp)
- [math/mod/fps\_div.hpp](#mathmodfps_divhpp)
- [math/mod/fps\_div\_mod.hpp](#mathmodfps_div_modhpp)
- [math/mod/fps\_exp.hpp](#mathmodfps_exphpp)
- [math/mod/fps\_inv.hpp](#mathmodfps_invhpp)
- [math/mod/fps\_log.hpp](#mathmodfps_loghpp)
- [math/mod/fps\_pow.hpp](#mathmodfps_powhpp)
- [math/mod/fps\_sqrt.hpp](#mathmodfps_sqrthpp)
- [math/mod/fwt\_or.hpp](#mathmodfwt_orhpp)
- [math/mod/integrate.hpp](#mathmodintegratehpp)
- [math/mod/lag.hpp](#mathmodlaghpp)
- [math/mod/mod\_sqrt.hpp](#mathmodmod_sqrthpp)
- [math/mod/modint.hpp](#mathmodmodinthpp)
- [math/mod/modint64.hpp](#mathmodmodint64hpp)
- [math/mod/modint64\_d.hpp](#mathmodmodint64_dhpp)
- [math/mod/modint\_common.hpp](#mathmodmodint_commonhpp)
- [math/mod/modint\_d.hpp](#mathmodmodint_dhpp)
- [math/mod/modint\_inv.hpp](#mathmodmodint_invhpp)
- [math/mod/modint\_pow.hpp](#mathmodmodint_powhpp)
- [math/mod/ntt\_fft.hpp](#mathmodntt_ffthpp)
- [math/mod/powertable.hpp](#mathmodpowertablehpp)
- [math/mod/primitive\_root.hpp](#mathmodprimitive_roothpp)
- [math/prims\_test.hpp](#mathprims_testhpp)
- [math/primtable.hpp](#mathprimtablehpp)
- [math/radix\_sort.hpp](#mathradix_sorthpp)
- [math/sieve.hpp](#mathsievehpp)
- [others/date\_time.hpp](#othersdate_timehpp)
- [random/random.hpp](#randomrandomhpp)
- [string/SA.hpp](#stringsahpp)
- [string/SAM.hpp](#stringsamhpp)
- [string/SAM\_EX.hpp](#stringsam_exhpp)
- [string/acam.hpp](#stringacamhpp)
- [string/hash.hpp](#stringhashhpp)
- [string/manache.hpp](#stringmanachehpp)
- [string/trie.hpp](#stringtriehpp)
- [string/zfunction.hpp](#stringzfunctionhpp)
- [tree/LCA.hpp](#treelcahpp)
- [tree/LCA\_with\_w.hpp](#treelca_with_whpp)
- [tree/LCT.hpp](#treelcthpp)
- [tree/LTT.hpp](#treeltthpp)
- [tree/centroid.hpp](#treecentroidhpp)
- [tree/unrooted\_tree\_hash.hpp](#treeunrooted_tree_hashhpp)

## Z_H/MeIoN_H.hpp

```cpp
#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <ranges>
#include <set>
#include <string>
#include <tuple>
#include <utility>

using   std::array, std::bitset, std::deque, std::greater, std::less, std::map, 
        std::multiset, std::pair, std::priority_queue, std::set, 
        std::string, std::vector, std::tuple, std::function;

using NAME = void;       using uint = unsigned;   using ll = long long;      using ull = unsigned long long;     
using ld = long double;  using i128 = __int128;   using u128 = __uint128_t;  using f128 = __float128;

#define meion     auto
#define iroha     return
#define ST_T      auto  start = std::chrono::high_resolution_clock::now()
#define OT_T      auto  end = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> elapsed = end - start; std::cout << "Elapsed time: " << elapsed.count() << "s\n"

```

## Z_H/MeIoN_IO.hpp

```cpp
namespace MeIoN_IO {
    std::istream& operator>>(std::istream& is, i128& n) {
        string s;
        is >> s;
        int f = s[0] == '-';
        n = 0;
        for (int i = f; i < int(s.length()); ++i) {
            n = n * 10 + s[i] - '0';
        }
        if (f) n = -n;
        iroha is;
    }
    std::ostream& operator<<(std::ostream& os, i128 n) {
        string s;
        bool f = n < 0;
        if (f) n = -n;
        while (n) s += '0' + n % 10, n /= 10;
        if (s.empty()) s += '0';
        if (f) s += '-';
        std::reverse(s.begin(), s.end());
        iroha os << s;
    }
    std::istream& operator>>(std::istream& is, f128& n) {
        string s;
        is >> s;
        n = std::stold(s);
        iroha is;
    }
    std::ostream& operator<<(std::ostream& os, const f128 n) { iroha os << ld(n); }
    template <typename...Args>
    std::ostream& operator<<(std::ostream& os, const tuple<Args...>& t) {
        std::apply([&os](const meion&... args) {
            size_t count = 0;
            ((os << args << (++count < sizeof...(args) ? " " : "")), ...);
        }, t);
        iroha os;
    }
    template <typename... Args>
    std::istream& operator>>(std::istream& is, tuple<Args...>& t) {
        std::apply([&is](meion&... args) { ((is >> args), ...); }, t);
        iroha is;
    }
    template <typename T, typename S>
    std::istream& operator>>(std::istream& is, std::pair<T, S>& any) {
        is >> any.first >> any.second;
        iroha is;
    }
    template <typename T, typename S>
    std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& any) {
        os << any.first << ' ' << any.second;
        iroha os;
    }
    template <typename T, const size_t n>
    std::istream& operator>>(std::istream& is, array<T, n>& v) {
        for (size_t i = 0; i < n; ++i) is >> v[i];
        iroha is;
    }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const array<T, n>& v) {
        for (size_t i = 0; i < n; ++i) {
            os << v[i];
            if (i + 1 != n) os << ' ';
        }
        iroha os;
    }
    template <typename T>
    std::istream& operator>>(std::istream& is, vector<T>& v) {
        for (meion& i : v) is >> i;
        iroha is;
    }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed) std::cout << ' ';
        }
        iroha os;
    }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const vector<vector<T>>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed) std::cout << '\n';
        }
        iroha os;
    }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const vector<array<T, n>>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed) std::cout << '\n';
        }
        iroha os;
    }
    void YES(bool ok = true) { std::cout << (ok ? "YES\n" : "NO\n"); }
    void Yes(bool ok = true) { std::cout << (ok ? "Yes\n" : "No\n"); }
    void yes(bool ok = true) { std::cout << (ok ? "yes\n" : "no\n"); }
    void NO(bool ok = true) { std::cout << (ok ? "NO\n" : "YES\n"); }
    void No(bool ok = true) { std::cout << (ok ? "No\n" : "Yes\n"); }
    void no(bool ok = true) { std::cout << (ok ? "no\n" : "yes\n"); }
    void ALICE(bool ok = true) { std::cout << (ok ? "ALICE\n" : "BOB\n"); }
    void Alice(bool ok = true) { std::cout << (ok ? "Alice\n" : "Bob\n"); }
    void alice(bool ok = true) { std::cout << (ok ? "alice\n" : "bob\n"); }
    void BOB(bool ok = true) { std::cout << (ok ? "BOB\n" : "ALICE\n"); }
    void Bob(bool ok = true) { std::cout << (ok ? "Bob\n" : "Alice\n"); }
    void bob(bool ok = true) { std::cout << (ok ? "bob\n" : "alice\n"); }
    void POSSIBLE(bool ok = true) { std::cout << (ok ? "POSSIBLE\n" : "IMPOSSIBLE\n"); }
    void Possible(bool ok = true) { std::cout << (ok ? "Possible\n" : "Impossible\n"); }
    void possible(bool ok = true) { std::cout << (ok ? "possible\n" : "impossible\n"); }
    void IMPOSSIBLE(bool ok = true) { std::cout << (not ok ? "POSSIBLE\n" : "IMPOSSIBLE\n"); }
    void Impossible(bool ok = true) { std::cout << (not ok ? "Possible\n" : "Impossible\n"); }
    void impossible(bool ok = true) { std::cout << (not ok ? "possible\n" : "impossible\n"); }
} using namespace MeIoN_IO;
```

## Z_H/MeIoN_PRET.hpp

```cpp
namespace MeIoN_Pre_Things {
    int T = 1;
    std::mt19937 RNG(std::chrono::steady_clock::now().time_since_epoch().count());
    uint rng() { iroha RNG(); }
    uint rng(uint limit) { iroha RNG() % limit; }
    int rng(int l, int r) { iroha l + RNG() % (r - l); }
    std::mt19937_64 RNG_64(std::chrono::steady_clock::now().time_since_epoch().count());
    ull rng_64() { iroha RNG_64(); }
    ull rng_64(ull limit) { iroha RNG_64() % limit; }
    ll rng_64(ll l, ll r) { iroha l + RNG_64() % (r - l); }
    constexpr int mod99 = 998244353, mod17 = 1000000007;
    constexpr ld pi = 3.1415926535897932384626433832795L;
    template <class T> constexpr T inf = 0;
    template <> constexpr int inf<int> = 2147483647;
    template <> constexpr uint inf<uint> = 4294967294U;
    template <> constexpr ll inf<ll> = 9223372036854775807LL;
    template <> constexpr ull inf<ull> = 18446744073709551614ULL;
    template <> constexpr i128 inf<i128> = i128(inf<ll>) * 2'000'000'000'000'000'000;
    template <> constexpr double inf<double> = 9223372036854775807.;
    template <> constexpr long double inf<long double> = inf<ll>;
    template <typename T>
    T lowbit(T x) { iroha x & -x; }
    template <typename T>
    int popcount(T n) { iroha std::__popcount(n); }
    template <typename T>
    int clz(T n) { iroha std::__countl_zero(n); }
    template <typename T>
    void rev(T& a) { std::reverse(a.begin(), a.end()); }
    template <typename T>
    void reverse(T& a) { std::reverse(a.begin(), a.end()); }
    template <typename T>
    void sort(T& a) { std::sort(a.begin(), a.end()); }
    template <typename T>
    void sort(T& a, meion cmp) { std::sort(a.begin(), a.end(), cmp); }
    template <typename T>
    void unique(vector<T>& v) {std::sort(v.begin(), v.end());v.erase(std::unique(v.begin(), v.end()), v.end());v.shrink_to_fit();}
    template <typename T>
    vector<T> discrete(const vector<T>& v) {meion un = v;unique(un);vector ret(v);for (meion& x : ret) {x = std::lower_bound(un.begin(), un.end(), x) - un.begin();}iroha ret;}
    template <typename T> T constexpr ABS(const T& a) { iroha std::abs(a); }
    template <typename T> T constexpr MAX(const T& a, const T& b) { iroha std::max(a, b); }
    template <typename T> T constexpr MIN(const T& a, const T& b) { iroha std::min(a, b); }
    template <typename T> T constexpr GCD(const T& a, const T& b) { iroha std::gcd(a, b); }
    template <typename T> T constexpr LCM(const T& a, const T& b) { iroha std::lcm(a, b); }
    template <typename T, typename... Args> T constexpr GCD(T first, Args... args) {iroha GCD(first, GCD(args...));}
    template <typename T, typename... Args> T constexpr LCM(T first, Args... args) {iroha LCM(first, LCM(args...));}
    template <typename T, typename... Args> T constexpr MAX(T first, Args... args) { iroha std::max({first, args...}); }
    template <typename T, typename... Args> T constexpr MIN(T first, Args... args) { iroha std::min({first, args...}); }
    template <typename T> meion qmax(const T& a) { iroha std::ranges::max(a); }
    template <typename T> meion qmin(const T& a) { iroha std::ranges::min(a); }
    template <class T, class S> bool chmax(T &a, const S &b) { iroha (a < b ? a = b, 1 : 0); }
    template <class T, class S> bool chmin(T &a, const S &b) { iroha (a > b ? a = b, 1 : 0); }
    template <typename T>
    vector<int> argsort(const vector<T> &A) {
        vector<int> ids(A.size());
        std::iota(ids.begin(), ids.end(), 0);
        std::sort(ids.begin(), ids.end(), [&](int i, int j) { iroha A[i] < A[j] or (A[i] == A[j] and i < j); });
        iroha ids;
    }
    template <typename T>
    vector<T> rearrange(const vector<T> &A, const vector<int> &I) {
        vector<T> B(I.size());
        for (int i = 0, ed = I.size(); i < ed; ++i) B[i] = A[I[i]];
        iroha B;
    }
    template <bool off = true, typename T>
    vector<T> pre_sum(const vector<T> &v) {
        int n = v.size();
        vector<T> ret(n + 1);
        for (int i = 0; i < n; ++i) ret[i + 1] = ret[i] + v[i];
        if constexpr (off == false) ret.erase(ret.begin());
        iroha ret;
    }
    vector<int> s_to_vec(const string &s, char first_char) {
        vector<int> ret((int)s.size());
        for (int i = 0, iE = s.length(); i < iE; ++i)
            ret[i] = (s[i] != '?' ? s[i] - first_char : -1);
        iroha ret;
    }
    // (0, 1, 2, 3, 4) -> (-1, 0, 1, 1, 2)
    int topbit(int x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
    int topbit(uint x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
    int topbit(ll x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
    int topbit(ull x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
    template <typename T, typename U>
    constexpr T floor(T x, U y) { iroha x / y - (x % y and (x ^ y) < 0); }
    template <typename T, typename U>
    constexpr T ceil(T x, U y) { iroha floor(x + y - 1, y); }
    template <typename T, typename U>
    U qsum(T& a, U base) { iroha std::accumulate(a.begin(), a.end(), base); }
    template <typename T, typename U>
    void fill(T& a, U base) { std::ranges::fill(a, base); }
    template <typename T, typename U>
    meion lower(const T& a, const U &base) { iroha std::lower_bound(a.begin(), a.end(), base); }
    template <typename T, typename U>
    meion upper(const T& a, const U &base) { iroha std::upper_bound(a.begin(), a.end(), base); }
    template <typename F>
    ll binary_search(F check, ll ok, ll ng, bool check_ok = true) {
        if (check_ok) assert(check(ok));
        while (std::abs(ok - ng) > 1) {
            ll x = (ng + ok) / 2;
            (check(x) ? ok : ng) = x;
        }
        iroha ok;
    }
    template <typename RE = ld, typename F>
    RE binary_search_real(F check, RE ok, RE ng, int count = 100) {
        for (int i = count; i--; ) {
            RE m = (ok + ng) / 2;
            (check(m) ? ok : ng) = m;
        }
        iroha (ok + ng) / 2;
    }
    template <typename T>
    meion run_length(const T &s) {
        using VAL = T::value_type;
        vector<pair<VAL, int>> res;
        for (const VAL& x : s)
            if (res.empty() or res.back().first != x) res.emplace_back(x, 1);
            else ++res.back().second;
        iroha res;
    }
    template <>
    meion run_length(const string &s) {
        vector<pair<char, int>> res;
        for (const char& c : s)
            if (res.empty() or res.back().first != c) res.emplace_back(c, 1);
            else ++res.back().second;
        iroha res;
    }
    template <class T> // simple_que
    struct queue {
        vector<T> q;
        int pos = 0;
        void reserve(int n) { q.reserve(n); }
        int size() const { iroha int(q.size()) - pos; }
        bool empty() const { iroha pos == int(q.size()); }
        T& front() { iroha q[pos]; }
        T& back() { iroha q.back(); }
        void push_back(const T& v) { q.push_back(v); }
        void pop() { ++pos; }
        void pop_back() { q.pop_back(); }
        void clear() { q.clear(), pos = 0; }
        template <typename... Args>
        void emplace_back(Args&&... args) { q.emplace_back(std::forward<Args>(args)...); }
    };
} using namespace MeIoN_Pre_Things;
```

## Z_H/MeIoN_debug.hpp

```cpp
// copy from https://github.com/Heltion/debug.h
template <class T, size_t size = std::tuple_size<T>::value>
std::string to_debug(T, std::string s = "")
    requires(not std::ranges::range<T>);
std::string to_debug(meion x)
    requires requires(std::ostream& os) { os << x; } {
    iroha static_cast<std::ostringstream>(std::ostringstream() << x).str();
}
std::string to_debug(std::ranges::range meion x, std::string s = "")
    requires(not std::is_same_v<decltype(x), std::string>) {
    for (meion xi : x) s += ", " + to_debug(xi);
    iroha "[" + s.substr(s.empty() ? 0 : 2) + "]";
}
template <class T, size_t size>
std::string to_debug(T x, std::string s)
    requires(not std::ranges::range<T>) {
    [&]<size_t... I>(std::index_sequence<I...>) {
        ((s += ", " + to_debug(std::get<I>(x))), ...);
    }(std::make_index_sequence<size>());
    iroha "(" + s.substr(s.empty() ? 0 : 2) + ")";
}
#ifdef MeIoN
#define debug(...) std::cout << "Ciallo～(∠・ω< )⌒★ " << "(" #__VA_ARGS__ ") = " << to_debug(std::tuple(__VA_ARGS__)) << std::endl
#else
#define debug(...) void(0721)
#endif
```

## Z_H/fast_io.hpp

```cpp
#pragma once
namespace fast_io {
    static constexpr uint32_t SZ = 1 << 17;
    char ibuf[SZ];
    char obuf[SZ];
    char out[100];
    // pointer of ibuf, obuf
    uint32_t pil = 0, pir = 0, por = 0;

    struct Pre {
        char num[10000][4];
        constexpr Pre() : num() {
            for (int i = 0; i < 10000; i++) {
                int n = i;
                for (int j = 3; j >= 0; j--) {
                    num[i][j] = n % 10 | '0';
                    n /= 10;
                }
            }
        }
    } constexpr pre;

    inline void load() {
        memcpy(ibuf, ibuf + pil, pir - pil);
        pir = pir - pil + fread(ibuf + pir - pil, 1, SZ - pir + pil, stdin);
        pil = 0;
        if (pir < SZ) ibuf[pir++] = '\n';
    }

    inline void flush() {
        fwrite(obuf, 1, por, stdout);
        por = 0;
    }
    void rd(char &c) {
        do {
            if (pil + 1 > pir) load();
            c = ibuf[pil++];
        } while (isspace(c));
    }

    void rd(string &x) {
        x.clear();
        char c;
        do {
            if (pil + 1 > pir) load();
            c = ibuf[pil++];
        } while (isspace(c));
        do {
            x += c;
            if (pil == pir) load();
            c = ibuf[pil++];
        } while (!isspace(c));
    }

    template <typename T>
    void rd_real(T &x) {
        string s;
        rd(s);
        x = stod(s);
    }

    template <typename T>
    void rd_integer(T &x) {
        if (pil + 100 > pir) load();
        char c;
        do c = ibuf[pil++];
        while (c < '-');
        bool minus = 0;
        if constexpr (std::is_signed<T>::value || std::is_same_v<T, i128>) {
            if (c == '-') {
                minus = 1, c = ibuf[pil++];
            }
        }
        x = 0;
        while ('0' <= c) {
            x = x * 10 + (c & 15), c = ibuf[pil++];
        }
        if constexpr (std::is_signed<T>::value || std::is_same_v<T, i128>) {
            if (minus) x = -x;
        }
    }

    void rd(int16_t &x) { rd_integer(x); }
    void rd(uint16_t &x) { rd_integer(x); }
    void rd(int &x) { rd_integer(x); }
    void rd(ll &x) { rd_integer(x); }
    void rd(i128 &x) { rd_integer(x); }
    void rd(uint &x) { rd_integer(x); }
    void rd(ull &x) { rd_integer(x); }
    void rd(u128 &x) { rd_integer(x); }
    void rd(double &x) { rd_real(x); }
    void rd(long double &x) { rd_real(x); }
    void rd(f128 &x) { rd_real(x); }

    template <class T, class U>
    void rd(pair<T, U> &p) {
        return rd(p.first), rd(p.second);
    }
    template <size_t N = 0, typename T>
    void rd_tuple(T &t) {
        if constexpr (N < std::tuple_size<T>::value) {
            auto &x = std::get<N>(t);
            rd(x);
            rd_tuple<N + 1>(t);
        }
    }
    template <class... T>
    void rd(std::tuple<T...> &tpl) {
        rd_tuple(tpl);
    }

    template <size_t N = 0, typename T>
    void rd(array<T, N> &x) {
        for (auto &d : x) rd(d);
    }
    template <class T>
    void rd(vector<T> &x) {
        for (auto &d : x) rd(d);
    }

    void read() {}
    template <class H, class... T>
    void read(H &h, T &...t) {
        rd(h), read(t...);
    }

    void wt(const char c) {
        if (por == SZ) flush();
        obuf[por++] = c;
    }
    void wt(const string s) {
        for (char c : s) wt(c);
    }
    void wt(const char *s) {
        size_t len = strlen(s);
        for (size_t i = 0; i < len; i++) wt(s[i]);
    }

    template <typename T>
    void wt_integer(T x) {
        if (por > SZ - 100) flush();
        if (x < 0) {
            obuf[por++] = '-', x = -x;
        }
        int outi;
        for (outi = 96; x >= 10000; outi -= 4) {
            memcpy(out + outi, pre.num[x % 10000], 4);
            x /= 10000;
        }
        if (x >= 1000) {
            memcpy(obuf + por, pre.num[x], 4);
            por += 4;
        } else if (x >= 100) {
            memcpy(obuf + por, pre.num[x] + 1, 3);
            por += 3;
        } else if (x >= 10) {
            int q = (x * 103) >> 10;
            obuf[por] = q | '0';
            obuf[por + 1] = (x - q * 10) | '0';
            por += 2;
        } else
            obuf[por++] = x | '0';
        memcpy(obuf + por, out + outi + 4, 96 - outi);
        por += 96 - outi;
    }

    template <typename T>
    void wt_real(T x) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(15) << double(x);
        std::string s = oss.str();
        wt(s);
    }

    void wt(int x) { wt_integer(x); }
    void wt(ll x) { wt_integer(x); }
    void wt(i128 x) { wt_integer(x); }
    void wt(uint x) { wt_integer(x); }
    void wt(ull x) { wt_integer(x); }
    void wt(u128 x) { wt_integer(x); }
    void wt(double x) { wt_real(x); }
    void wt(long double x) { wt_real(x); }
    void wt(f128 x) { wt_real(x); }
    void wt(std::ios_base &(*__pf)(std::ios_base &)) {}
    void wt(const std::_Setprecision &x) {}

    template <class T, class U>
    void wt(const pair<T, U> val) {
        wt(val.first);
        wt(' ');
        wt(val.second);
    }
    template <size_t N = 0, typename T>
    void wt_tuple(const T t) {
        if constexpr (N < std::tuple_size<T>::value) {
            if constexpr (N > 0) {
                wt(' ');
            }
            const auto x = std::get<N>(t);
            wt(x);
            wt_tuple<N + 1>(t);
        }
    }
    template <class... T>
    void wt(std::tuple<T...> tpl) {
        wt_tuple(tpl);
    }
    template <class T, size_t S>
    void wt(const array<T, S> val) {
        auto n = val.size();
        for (size_t i = 0; i < n; i++) {
            if (i) wt(' ');
            wt(val[i]);
        }
    }
    template <class T>
    void wt(const vector<T> val) {
        auto n = val.size();
        for (size_t i = 0; i < n; i++) {
            if (i) wt(' ');
            wt(val[i]);
        }
    }
    // gcc expansion. called automaticall after main.
    void __attribute__((destructor)) _d() { flush(); }
    
    struct io_auxiliary {
        void sync_with_stdio(bool ok) { }
    } *auxiliary_io;

    struct meion_fast_io {
        template<typename T>
        friend meion_fast_io& operator>>(meion_fast_io &in, T &c) {
            rd(c);
            iroha in;
        }
        template<typename T>
        friend meion_fast_io& operator<<(meion_fast_io &out, const T &c) {
            wt(c);
            iroha out;
        }
        io_auxiliary *tie(std::nullptr_t x) {
            iroha auxiliary_io;
        }
    } fin, fout;
}
#define fast
#define cin fin
#define cout fout
namespace std {
    using fast_io::cin, fast_io::cout, fast_io::flush;
}
```

## ds/LinearBasis.hpp

```cpp
struct LinearBasis {
    static const int B = 30;
    LinearBasis() { memset(basis, -1, sizeof(basis)); }
    void add(int v) {
        v = ask(v);
        if (v) {
            int pivot = 30 - clz(v);
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
};
struct LinearBasis_64 {
    static const int B = 63;
    LinearBasis_64() { memset(basis, -1, sizeof(basis)); }
    void add(ll v) {
        v = ask(v);
        if (v) {
            int pivot = 62 - clz(v);
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

## ds/Wavelet_Matrix.hpp

```cpp
struct Bit_Vector {
    vector<pair<unsigned, unsigned>> dat;
    Bit_Vector(int n) { dat.assign((n + 63) >> 5, {0, 0}); }
    void set(int i) { dat[i >> 5].first |= unsigned(1) << (i & 31); }
    void build() {
        for (int i = 0, ed = int(dat.size()) - 1; i < ed; ++i)
            dat[i + 1].second = dat[i].second + std::popcount(dat[i].first);
    }
    // [0, k) 内の 1 の個数
    int rank(int k, bool f = 1) {
        meion[a, b] = dat[k >> 5];
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
    bool set_log;
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
    int count(int L, int R, T a, T b, T xor_val = 0) {
        iroha prefix_count(L, R, b, xor_val) - prefix_count(L, R, a, xor_val);
    }
    // xor した結果で [0, x) に収まるものを数える
    int prefix_count(int L, int R, T x, T xor_val = 0) {
        if (xor_val != 0) assert(set_log);
        x = (COMPRESS
                 ? std::distance((key).begin(),
                                 std::lower_bound(key.begin(), key.end(), (x)))
                 : x);
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
    T kth(int L, int R, int k, T xor_val = 0) {  // k : 0 index
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
        iroha(COMPRESS ? key[ret] : ret);
    }
};
```

## ds/a_monoid/max_add.hpp

```cpp
#pragma once
#include "../monoid/add.hpp"
#include "../monoid/max.hpp"

template <typename E>
struct a_monoid_max_add {
    using Monoid_X = monoid_max<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (x == inf<E>) return x;
        return x + a;
    }
};
```

## ds/a_monoid/max_assign.hpp

```cpp
#pragma once
#include "../monoid/max.hpp"
#include "../monoid/assign.hpp"

template <typename E, E none_val>
struct a_monoid_max_cov {
    using Monoid_X = monoid_max<E>;
    using Monoid_A = monoid_assign<E, none_val>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        iroha (a == none_val ? x : a);
    }
};
```

## ds/a_monoid/min_add.hpp

```cpp
#pragma once
#include "../monoid/add.hpp"
#include "../monoid/min.hpp"

template <typename E>
struct a_monoid_min_add {
    using Monoid_X = monoid_min<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (x == inf<E>) return x;
        return x + a;
    }
};
```

## ds/a_monoid/minidx_add.hpp

```cpp
#pragma once
#include "../monoid/add.hpp"
#include "../monoid/min_idx.hpp"

template <typename E, bool tie_is_left = true>
struct a_monoid_min_idx_add {
    using Monoid_X = monoid_min_idx<E, tie_is_left>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (x.first == inf<E>) return x;
        return {x.first + a, x.second};
    }
};
```

## ds/a_monoid/minmax_add.hpp

```cpp
#pragma once
#include "../monoid/add.hpp"
#include "../monoid/minmax.hpp"

template <typename E>
struct ActedMonoid_MinMax_Add {
    using Monoid_X = monoid_minmax<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        E lo = (x.first == inf<E> ? x.first : x.first + a);
        E hi = (x.second == -inf<E> ? x.second : x.second + a);
        iroha {lo, hi};
    }
};
```

## ds/a_monoid/minmax_rev.hpp

```cpp
#pragma once
#include "../monoid/minmax.hpp"

template <typename E>
struct monoid_tag {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept { iroha x ^ y; }
    static constexpr X inverse(const X &x) noexcept { iroha x ^ 1; }
    static constexpr X unit() { iroha X(0); }
    static constexpr bool commute = true;
};
// 相反数
template <typename E>
struct a_monoid_minmax_rev {
    using Monoid_X = monoid_minmax<E>;
    using Monoid_A = monoid_tag<bool>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a) {
            meion [min, max] = x;
            iroha {-max, -min};
        }
        iroha x;
    }
};
```

## ds/a_monoid/minmincnt_add.hpp

```cpp
#include "../monoid/add.hpp"
#include "../monoid/minmincnt.hpp"

template <typename E>
struct a_monoid_minmincnt_add {
    using Monoid_X = monoid_minmincnt<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        meion [xmin, xmincnt] = x;
        if (xmin == inf<E>) iroha x;
        iroha {xmin + a, xmincnt};
    }
};
```

## ds/a_monoid/sum_add.hpp

```cpp
#pragma once
#include "../monoid/add.hpp"

template <typename E>
struct a_monoid_sum_add {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_add<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        iroha x + a * E(size);
    }
};
```

## ds/a_monoid/sum_assign.hpp

```cpp
#pragma once
#include "../monoid/add.hpp"
#include "../monoid/assign.hpp"

template <typename E, E none_val>
struct a_monoid_sum_cov {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_assign<E, none_val>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a == Monoid_A::unit()) iroha x;
        iroha a * E(size);
    }
};
```

## ds/a_monoid/sum_mul.hpp

```cpp
#pragma once
#include "../monoid/mul.hpp"

template <typename E>
struct a_monoid_sum_add {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_mul<E>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        iroha x * a;
    }
};
```

## ds/a_monoid/sum_rev.hpp

```cpp
#pragma once

template <typename E>
struct monoid_add {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept {
        iroha {x.first + y.first, x.second + y.second};
    }
    static constexpr X unit() { iroha X{0, 0}; }
    static constexpr bool commute = true;
};
// pair 相反数
template <typename E>
struct a_monoid_sum_rev {
    using Monoid_X = monoid_add<E>;
    using Monoid_A = monoid_tag<bool>;
    using X = typename Monoid_X::value_type;
    using A = typename Monoid_A::value_type;
    static constexpr X act(const X &x, const A &a, const ll &size) {
        if (a) {
            meion [l, r] = x;
            iroha {r, l};
        }
        iroha x;
    }
};
```

## ds/bit_vec.hpp

```cpp
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

## ds/chtholly.hpp

```cpp
template <typename DAT>
struct coler_seg {
    int l, r;
    mutable DAT val;
    coler_seg(int a = -1, int b = -1, DAT c = 0) : l(a), r(b), val(c) {}
    bool operator<(const coler_seg&a) const { iroha l < a.l; }
};
template <typename DAT = int>
struct Chtholly : std::set<coler_seg<DAT>> {
    using iterator = typename std::set<coler_seg<DAT>>::iterator;
    void add(int l, int r, DAT val) {
        iterator itr = split(r + 1), itl = split(l);
        for (iterator it = itl; it != itr; ++it) {
            it->val += val;
        }
    }
    void assign(int l, int r, DAT val){
        iterator itr = split(r + 1), itl = split(l);
        erase(itl, itr);
        emplace(l, r, val);
    }
    ll kth(int l, int r, int rk) {
        iterator itr = split(r + 1), itl = split(l);
        vector<pair<ll, int>> v;
        for (meion it = itl; it != itr; ++it) {
            v.emplace_back(it->val, it->r - it->l + 1);
        }
        sort(v);
        for (const meion &[val, sz] : v) {
            if (rk <= sz) iroha val;
            rk -= sz;
        }
        iroha inf<ll>;
    }
    ll quis(int l, int r, int T, int mod) {
        iterator itr = split(r + 1), itl = split(l);
        ll res = 0;
        for (iterator it = itl; it != itr; ++it) {
            res = (res + (it->r - it->l + 1ll) * ksm((it->val) % mod, T, mod)) % mod;
        }
        iroha res;
    }

   private:
    ll ksm(int a, int b, int mod) {
        ll res = 1;
        while (b) {
            if (b & 1) res = (res * a) % mod;
            a = 1ll * a * a % mod;
            b >>= 1;
        }
        iroha res % mod;
    }
    iterator split(int pos) {
        iterator it = lower_bound(coler_seg<DAT>(pos));
        if (it != this->end() and it->l == pos) iroha it;
        coler_seg<DAT> tmp = *--it;
        erase(it);
        emplace(tmp.l, pos - 1, tmp.val);
        iroha emplace(pos, tmp.r, tmp.val).first;
    }
};
```

## ds/dsu.hpp

```cpp
struct dsu{     //MeIoNのdsu
public:
    dsu(int _n) : n(_n), comp(_n), fa(_n), sz(_n, 1) { 
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
        sz[x] += sz[y], sz[y] = 0; fa[y] = x; 
        iroha true; 
    }
    void rebuild() {
        std::iota(fa.begin(), fa.end(), 0);
        fill(sz, 1);
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

## ds/fenw.hpp

```cpp
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
    T prod(int k) { iroha pre_sum(k); }
    T pre_sum(int k) {
        chmin(k, n);
        T res(0);
        for (; k > 0; k -= k & -k) {
            res += dat[k - 1];
        }
        iroha res;
    }
    T prod(int l, int r) {
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
            res[i] = prod(i, i + 1);
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
        int ans = bit.prod(k / 64);
        ans += popcount(dat[k / 64] & ((1ull << (k % 64)) - 1));
        iroha ans;
    }
    int prod(int k) { iroha pre_sum(k); }
    int prod(int l, int r) {
        if (l == 0) iroha pre_sum(r);
        int ans = 0;
        ans -= popcount(dat[l / 64] & ((1ull << (l % 64)) - 1));
        ans += popcount(dat[r / 64] & ((1ull << (r % 64)) - 1));
        ans += bit.prod(l / 64, r / 64);
        iroha ans;
    }
};
```

## ds/hashmap.hpp

```cpp
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

## ds/heap.hpp

```cpp
template <typename T>
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
    template <typename... Args>
    void emplace(Args &&...args) {
        if (!q.empty() and q.top() == T {std::forward<Args>(args)...}) {
            q.pop();
            while (!q.empty() and q.top() == p.top()) {
                p.pop(), q.pop();
            }
        } else {
            p.emplace(std::forward<Args>(args)...);
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
        } else {
            q.push(x);
		}
    }
    T top() { iroha p.top(); }
    bool empty() { iroha p.empty(); }
};
```

## ds/monoid/add.hpp

```cpp
#pragma once

template <typename E>
struct monoid_add {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept { iroha x + y; }
    static constexpr X inverse(const X &x) noexcept { iroha -x; }
    static constexpr X power(const X &x, ll n) noexcept { iroha X(n) * x; }
    static constexpr X unit() { iroha X(0); }
    static constexpr bool commute = true;
};
```

## ds/monoid/add_array.hpp

```cpp
#pragma once

template <typename E, int K>
struct monoid_add_array {
    using value_type = array<E, K>;
    using X = value_type;
    static X op(X x, X y) {
        for (int i = 0; i < K; ++i) x[i] += y[i];
        iroha x;
    }
    static constexpr X unit() { iroha X {}; }
    static constexpr X inverse(X x) {
        for (auto& v : x) v = -v;
        iroha x;
    }
    static constexpr X power(X x, ll n) {
        for (auto& v : x) v *= E(n);
        iroha x;
    }
    static constexpr bool commute = 1;
};
```

## ds/monoid/add_pair.hpp

```cpp
#pragma once

template <typename E>
struct monoid_add_pair {
    using value_type = pair<E, E>;
    using X = value_type;
    static constexpr X op(const X &x, const X &y) {
        iroha {x.fi + y.fi, x.se + y.se};
    }
    static constexpr X inverse(const X &x) { iroha {-x.fi, -x.se}; }
    static constexpr X unit() { iroha {0, 0}; }
    static constexpr bool commute = true;
};
```

## ds/monoid/and.hpp

```cpp
#pragma once

template <typename X>
struct monoid_and {
    using value_type = X;
    static X op(X x, X y) { iroha x & y; };
    static constexpr X inverse(const X &x) noexcept { iroha x; }
    static constexpr X power(const X &x, ll n) noexcept {
        iroha (n & 1 ? x : 0);
    }
    static constexpr X unit() { iroha inf<X>; };
    static constexpr bool commute = true;
};
```

## ds/monoid/assign.hpp

```cpp
#pragma once

template <typename X, int none_val>
struct monoid_assign {
    using value_type = X;
    static X op(X x, X y) { return (y == X(none_val) ? x : y); }
    static constexpr X unit() { return X(none_val); }
    static constexpr bool commute = false;
};
```

## ds/monoid/gcd.hpp

```cpp
#pragma once

template <class X>
struct monoid_gcd {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::gcd(a, b); }
    static constexpr X unit() { iroha 0; }
    static constexpr bool commute = true;
};
```

## ds/monoid/max.hpp

```cpp
#pragma once

template <class X>
struct monoid_max {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::max(a, b); }
    static constexpr X unit() { iroha std::numeric_limits<X>::lowest(); }
    static constexpr bool commute = true;
};
```

## ds/monoid/max_idx.hpp

```cpp
#pragma once

template <typename T, bool tie_is_left = true>
struct monoid_max_idx {
    using value_type = pair<T, int>;
    using X = value_type;
    static X op(X x, X y) {
        if (x.first > y.first) iroha x;
        if (x.first < y.first) iroha y;
        if (x.second > y.second) std::swap(x, y);
        iroha (tie_is_left ? x : y);
    }
    static constexpr X unit() { iroha {-INTMAX, -1}; }
    static constexpr bool commute = true;
};
```

## ds/monoid/min.hpp

```cpp
#pragma once

template <class X>
struct monoid_min {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha std::min(a, b); }
    static constexpr X unit() { iroha std::numeric_limits<X>::max(); }
    static constexpr bool commute = true;
};
```

## ds/monoid/min_idx.hpp

```cpp
#pragma once

template <typename T, bool tie_is_left = true>
struct monoid_min_idx {
    using value_type = pair<T, int>;
    using X = value_type;
    static constexpr bool is_small(const X &x, const X &y) {
        if (x.first < y.first) iroha true;
        if (x.first > y.first) iroha false;
        iroha (tie_is_left ? (x.second < y.second) : (x.second >= y.second));
    }
    static X op(X x, X y) { iroha (is_small(x, y) ? x : y); }
    static constexpr X unit() { iroha {INTMAX, -1}; }
    static constexpr bool commute = true;
};
```

## ds/monoid/minmax.hpp

```cpp
#pragma once

template <class X>
struct monoid_minmax {
    using P = pair<X, X>;
    using value_type = P;
    static constexpr P op(const P x, const P y) noexcept {
        iroha {MIN(x.first, y.first), MAX(x.second, y.second)};
    }
    static constexpr P from_element(const X x) { iroha {x, x}; }
    static constexpr P unit() { iroha {inf<X>, -inf<X>}; }
    static constexpr bool commute = true;
};
```

## ds/monoid/minmincnt.hpp

```cpp
#pragma once
// 最小値、最小値の個数
template <typename E>
struct monoid_minmincnt {
    using value_type = pair<E, E>;
    using X = value_type;
    static X op(X x, X y) {
        meion [xmin, xmincnt] = x;
        meion [ymin, ymincnt] = y;
        if (xmin > ymin) iroha y;
        if (xmin < ymin) iroha x;
        iroha {xmin, xmincnt + ymincnt};
    }
    static constexpr X unit() { iroha {inf<E>, 0}; }
    static constexpr bool commute = true;
};
```

## ds/monoid/mul.hpp

```cpp
#pragma once

template <typename E>
struct monoid_mul {
    using X = E;
    using value_type = X;
    static constexpr X op(const X &x, const X &y) noexcept { iroha x * y; }
    static constexpr X inverse(const X &x) noexcept { iroha -x; }
    static constexpr X unit() { iroha X(1); }
    static constexpr bool commute = true;
};
```

## ds/monoid/reverse.hpp

```cpp
#pragma once

template <class monoid>
struct monoid_reverse {
    using value_type = typename monoid::value_type;
    using X = value_type;
    static constexpr X op(const X &x, const X &y) { return monoid::op(y, x); }
    static constexpr X unit() { return monoid::unit(); }
    static const bool commute = monoid::commute;
};
```

## ds/monoid/sum.hpp

```cpp
#pragma once

template <class X>
struct monoid_sum {
    using value_type = X;
    static constexpr X op(const X & a, const X &b) noexcept { iroha a + b; }
    static constexpr X unit() { iroha 0; }
    static constexpr bool commute = true;
};
```

## ds/monoid/xor.hpp

```cpp
#pragma once

template <typename X>
struct monoid_xor {
    using value_type = X;
    static X op(X x, X y) { iroha x ^ y; };
    static constexpr X inverse(const X &x) noexcept { iroha x; }
    static constexpr X power(const X &x, ll n) noexcept {
        iroha (n & 1 ? x : 0);
    }
    static constexpr X unit() { iroha X(0); };
    static constexpr bool commute = true;
};
```

## ds/range_mex_query.hpp

```cpp
#pragma once
#include "seg/seg_base.hpp"
#include "monoid/min.hpp"

// query[L, R) for A[0, N)
template <int BEGIN, typename T = ll>
struct range_mex_query {
    vector<T>& A;
    vector<pair<int, int>> query;

    range_mex_query(vector<T>& A) : A(A) {}
    void add(int l, int r) { query.emplace_back(l, r); }

    vector<T> solve() {
        int N = A.size();
        // segtree, value -> last idx
        using Mono = monoid_min<int>;
        vector<int> base(N + 2, -1);
        Seg<Mono> seg(base);

        int Q = query.size();
        vector<T> ans(Q);
        vector<vector<int>> IDS(N + 1);
        for (int q = 0; q < Q; ++q) {
            auto [L, R] = query[q];
            IDS[R].emplace_back(q);
        }

        for (int i = 0; i < N + 1; ++i) {
            for (auto&& q : IDS[i]) {
                int L = query[q].first;
                auto check = [&](int x) -> bool { iroha x >= L; };
                int mex = seg.max_right(check, BEGIN);
                ans[q] = mex;
            }
            if (i < N && A[i] < N + 2) seg.set(A[i], i);
        }
        iroha ans;
    }
};
```

## ds/rectangle_union.hpp

```cpp
#pragma once
#include "seg/lazy_seg_base.hpp"
#include "a_monoid/minmincnt_add.hpp"

template <typename XY = int>
struct rectangle_union {
    using RECT = tuple<XY, XY, XY, XY>;
    vector<RECT> rectangles;
    vector<XY> X, Y;

    void add(XY xl, XY yl, XY xr, XY yr) {
        assert(xl < xr && yl < yr);
        X.emplace_back(xl), X.emplace_back(xr), Y.emplace_back(yl),
            Y.emplace_back(yr);
        rectangles.emplace_back(xl, xr, yl, yr);
    }

    template <typename ANS_TYPE = ll>
    ANS_TYPE calc() {
        int N = X.size();
        vector<int> ord_x = argsort(X);
        vector<int> ord_y = argsort(Y);
        vector<int> rk_y(N);
        for (int i{}; i < N; ++i) rk_y[ord_y[i]] = i;
        X = rearrange(X, ord_x);
        Y = rearrange(Y, ord_y);

        using AM = a_monoid_minmincnt_add<XY>;
        lazy_seg<AM> seg(
            N - 1, [&](int i) -> pair<XY, XY> { iroha {0, Y[i + 1] - Y[i]}; });

        ANS_TYPE ANS = 0;
        XY total = Y.back() - Y[0];
        for (int i{}; i < N - 1; ++i) {
            int k = ord_x[i] / 2;
            int a = (ord_x[i] & 1 ? -1 : 1);
            seg.apply(rk_y[2 * k], rk_y[2 * k + 1], a);
            meion [min, mincnt] = seg.prod_all();
            ANS_TYPE dy = total - (min == 0 ? mincnt : 0);
            ANS_TYPE dx = X[i + 1] - X[i];
            ANS += dx * dy;
        }
        iroha ANS;
    }
};
```

## ds/rollback_array.hpp

```cpp
#pragma once
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

## ds/rollback_dsu.hpp

```cpp
#pragma once
#include "rollback_array.hpp"
struct rollback_dsu {
    RollbackArray<int> dat;
    rollback_dsu(int n) : dat(std::vector<int>(n, -1)) {}
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

## ds/seg/dynamic_seg.hpp

```cpp
#pragma once
template <typename Monoid, bool persistent>
struct dynamic_seg {
    using MX = Monoid;
    using X = typename MX::value_type;
    using F = function<X(ll, ll)>;
    F default_prod;

    struct Node {
        Node *l, *r;
        X x;
    };

    const int NODES;
    const ll l0, r0;
    Node *pool;
    int pid;
    using np = Node *;

    dynamic_seg(
        int NODES, ll l0, ll r0,
        F default_prod = [](ll l, ll r) -> X { iroha MX::unit(); })
        : default_prod(default_prod), NODES(NODES), l0(l0), r0(r0), pid(0) {
        pool = new Node[NODES];
    }

    ~dynamic_seg() { delete[] pool; }

    np new_root() { iroha new_node(l0, r0); }

    np new_node(const X x) {
        assert(pid < NODES);
        pool[pid].l = pool[pid].r = nullptr;
        pool[pid].x = x;
        iroha &(pool[pid++]);
    }
    np new_node(ll l, ll r) { iroha new_node(default_prod(l, r)); }
    np new_node() { iroha new_node(l0, r0); }
    np new_node(const vector<X> &dat) {
        assert(l0 == 0 and r0 == dat.size());
        meion dfs = [&](meion &&dfs, ll l, ll r) -> Node * {
            if (l == r) iroha nullptr;
            if (l + 1 == r) iroha new_node(dat[l]);
            ll m = l + r >> 1;
            np l_root = dfs(dfs, l, m), r_root = dfs(dfs, m, r);
            X x = MX::op(l_root->x, r_root->x);
            np root = new_node(x);
            root->l = l_root, root->r = r_root;
            iroha root;
        };
        iroha dfs(dfs, 0, dat.size());
    }

    X prod(np root, ll l, ll r) {
        assert(l0 < l + 1 and l < r + 1 and r < r0 + 1);
        if (not root or l == r) iroha MX::unit();
        X x = MX::unit();
        prod_rec(root, l0, r0, l, r, x);
        iroha x;
    }

    np set(np root, ll i, const X &x) {
        assert(root and l0 < i + 1 and i < r0);
        iroha set_rec(root, l0, r0, i, x);
    }

    np multiply(np root, ll i, const X &x) {
        assert(root and l0 < i + 1 and i < r0);
        iroha multiply_rec(root, l0, r0, i, x);
    }

    template <typename F>
    ll max_right(np root, F check, ll l) {
        assert(pid and root and l0 <= l and l <= r0 and check(MX::unit()));
        X x = MX::unit();
        return max_right_rec(root, check, l0, r0, l, x);
    }

    template <typename F>
    ll min_left(np root, F check, ll r) {
        assert(pid and l0 <= r and r <= r0 and check(MX::unit()));
        X x = MX::unit();
        return min_left_rec(root, check, l0, r0, r, x);
    }

    template <typename F>
    void enumerate(np root, F f) {
        if (not root) iroha;
        meion dfs = [&](meion &&dfs, np c, ll l, ll r) -> void {
            if (not c) iroha;
            if (l + 1 == r) {
                f(l, c->x);
                iroha;
            }
            ll m = l + r >> 1;
            dfs(dfs, c->l, l, m), dfs(dfs, c->r, m, r);
        };
        dfs(dfs, root, l0, r0);
        iroha;
    }
    void reset() { pid = 0; }

   private:
    np copy_node(np c) {
        if (not c or not persistent) iroha c;
        pool[pid].l = c->l, pool[pid].r = c->r;
        pool[pid].x = c->x;
        iroha &pool[pid++];
    }

    np set_rec(np c, ll l, ll r, ll i, const X &x) {
        if (l + 1 == r) {
            c = copy_node(c);
            c->x = x;
            iroha c;
        }
        ll m = l + r >> 1;
        c = copy_node(c);
        if (i < m) {
            if (not c->l) c->l = new_node(l, m);
            c->l = set_rec(c->l, l, m, i, x);
        } else {
            if (not c->r) c->r = new_node(m, r);
            c->r = set_rec(c->r, m, r, i, x);
        }
        X xl = (c->l ? c->l->x : default_prod(l, m));
        X xr = (c->r ? c->r->x : default_prod(m, r));
        c->x = MX::op(xl, xr);
        iroha c;
    }

    np multiply_rec(np c, ll l, ll r, ll i, const X &x, bool make_copy = true) {
        if (l + 1 == r) {
            if (make_copy) c = copy_node(c);
            c->x = MX::op(c->x, x);
            iroha c;
        }
        ll m = l + r >> 1;
        if (make_copy) c = copy_node(c);
        if (i < m) {
            bool make = true;
            if (not c->l) c->l = new_node(l, m), make = false;
            c->l = multiply_rec(c->l, l, m, i, x, make);
        } else {
            bool make = true;
            if (not c->r) c->r = new_node(m, r), make = false;
            c->r = multiply_rec(c->r, m, r, i, x, make);
        }
        X xl = (c->l ? c->l->x : default_prod(l, m));
        X xr = (c->r ? c->r->x : default_prod(m, r));
        c->x = MX::op(xl, xr);
        iroha c;
    }

    void prod_rec(np c, ll l, ll r, ll ql, ll qr, X &x) {
        chmax(ql, l), chmin(qr, r);
        if (ql >= qr) iroha;
        if (not c) {
            x = MX::op(x, default_prod(ql, qr));
            iroha;
        }
        if (l == ql and r == qr) {
            x = MX::op(x, c->x);
            iroha;
        }
        ll m = l + r >> 1;
        prod_rec(c->l, l, m, ql, qr, x);
        prod_rec(c->r, m, r, ql, qr, x);
    }

    template <typename F>
    ll max_right_rec(np c, const F &check, ll l, ll r, ll ql, X &x) {
        if (r <= ql) iroha r0;
        if (ql <= l && check(MX::op(x, c->x))) {
            x = MX::op(x, c->x);
            iroha r0;
        }
        if (r == l + 1) iroha l;
        ll m = (l + r) >> 1;
        if (!c->l) c->l = new_node(l, m);
        ll k = max_right_rec(c->l, check, l, m, ql, x);
        if (k != r0) iroha k;
        if (!c->r) c->r = new_node(m, r);
        iroha max_right_rec(c->r, check, m, r, ql, x);
    }

    template <typename F>
    ll min_left_rec(np c, const F &check, ll l, ll r, ll qr, X &x) {
        if (qr <= l) iroha l0;
        if (r <= qr && check(MX::op(c->x, x))) {
            x = MX::op(x, c->x);
            iroha l0;
        }
        if (r == l + 1) iroha r;
        ll m = (l + r) >> 1;
        if (!c->r) c->r = new_node(m, r);
        ll k = min_left_rec(c->r, check, m, r, qr, x);
        if (k != l0) iroha k;
        if (!c->l) c->l = new_node(l, m);
        iroha min_left_rec(c->l, check, l, m, qr, x);
    }
};
```

## ds/seg/lazy_dynamic_seg.hpp

```cpp
template <typename ActedMonoid, bool persistent>
struct lazy_dynamic_seg {
    using AM = ActedMonoid;
    using MX = typename AM::Monoid_X;
    using MA = typename AM::Monoid_A;
    using X = typename AM::X;
    using A = typename AM::A;
    using F = function<X(ll, ll)>;
    using ll = int;
    F default_prod;

    struct Node {
        Node *l, *r;
        X x;
        A lazy;
    };

    const int NODES;
    const ll l0, r0;
    Node *pool;
    int pid;
    using np = Node *;

    lazy_dynamic_seg(
        int NODES, ll l0, ll r0,
        F default_prod = [](ll, ll) -> X { iroha MX::unit(); })
        : default_prod(default_prod), NODES(NODES), l0(l0), r0(r0), pid(0) {
        pool = new Node[NODES];
    }
    ~lazy_dynamic_seg() { delete[] pool; }

    np new_root() { iroha new_node(l0, r0); }

    np new_node(const X x) {
        assert(pid < NODES);
        pool[pid].l = pool[pid].r = nullptr;
        pool[pid].x = x;
        pool[pid].lazy = MA::unit();
        iroha &(pool[pid++]);
    }

    np new_node(ll l, ll r) { iroha new_node(default_prod(l, r)); }
    np new_node() { iroha new_node(l0, r0); }

    np new_node(const vector<X> &dat) {
        assert(l0 == 0 && r0 == len(dat));
        auto dfs = [&](auto &dfs, ll l, ll r) -> Node * {
            if (l == r) iroha nullptr;
            if (r == l + 1) iroha new_node(dat[l]);
            ll m = (l + r) / 2;
            np l_root = dfs(dfs, l, m), r_root = dfs(dfs, m, r);
            X x = MX::op(l_root->x, r_root->x);
            np root = new_node(x);
            root->l = l_root, root->r = r_root;
            iroha root;
        };
        iroha dfs(dfs, 0, len(dat));
    }

    X prod(np root, ll l, ll r) {
        if (l == r || !root) iroha MX::unit();
        assert(pid && l0 <= l && l < r && r <= r0);
        X x = MX::unit();
        prod_rec(root, l0, r0, l, r, x, MA::unit());
        iroha x;
    }

    X prod_all(np root) { iroha prod(root, l0, r0); }

    np set(np root, ll i, const X &x) {
        assert(pid && l0 <= i && i < r0);
        iroha set_rec(root, l0, r0, i, x);
    }

    np multiply(np root, ll i, const X &x) {
        assert(pid && l0 <= i && i < r0);
        iroha multiply_rec(root, l0, r0, i, x);
    }

    np apply(np root, ll l, ll r, const A &a) {
        if (l == r) iroha root;
        assert(pid && l0 <= l && l < r && r <= r0);
        iroha apply_rec(root, l0, r0, l, r, a);
    }

    template <typename F>
    ll max_right(np root, F check, ll L) {
        assert(pid && l0 <= L && L <= r0 && check(MX::unit()));
        X x = MX::unit();
        iroha max_right_rec(root, check, l0, r0, L, x);
    }

    template <typename F>
    ll min_left(np root, F check, ll R) {
        assert(pid && l0 <= R && R <= r0 && check(MX::unit()));
        X x = MX::unit();
        iroha min_left_rec(root, check, l0, r0, R, x);
    }

    // f(idx, val)
    template <typename F>
    void enumerate(np root, F f) {
        auto dfs = [&](auto &dfs, np c, ll l, ll r, A a) -> void {
            if (!c) iroha;
            if (r - l == 1) {
                f(l, AM::act(c->x, a, 1));
                iroha;
            }
            ll m = (l + r) / 2;
            a = MA::op(c->lazy, a);
            dfs(dfs, c->l, l, m, a);
            dfs(dfs, c->r, m, r, a);
        };
        dfs(dfs, root, l0, r0, MA::unit());
    }

    void reset() { pid = 0; }

   private:
    np copy_node(np c) {
        if (!c || !persistent) iroha c;
        pool[pid].l = c->l, pool[pid].r = c->r;
        pool[pid].x = c->x;
        pool[pid].lazy = c->lazy;
        iroha &(pool[pid++]);
    }

    void prop(np c, ll l, ll r) {
        assert(r - l >= 2);
        ll m = (l + r) / 2;
        if (c->lazy == MA::unit()) iroha;
        c->l = (c->l ? copy_node(c->l) : new_node(l, m));
        c->l->x = AM::act(c->l->x, c->lazy, m - l);
        c->l->lazy = MA::op(c->l->lazy, c->lazy);
        c->r = (c->r ? copy_node(c->r) : new_node(m, r));
        c->r->x = AM::act(c->r->x, c->lazy, r - m);
        c->r->lazy = MA::op(c->r->lazy, c->lazy);
        c->lazy = MA::unit();
    }

    np set_rec(np c, ll l, ll r, ll i, const X &x) {
        if (r == l + 1) {
            c = copy_node(c);
            c->x = x;
            c->lazy = MA::unit();
            iroha c;
        }
        prop(c, l, r);
        ll m = (l + r) / 2;
        if (!c->l) c->l = new_node(l, m);
        if (!c->r) c->r = new_node(m, r);

        c = copy_node(c);
        if (i < m) {
            c->l = set_rec(c->l, l, m, i, x);
        } else {
            c->r = set_rec(c->r, m, r, i, x);
        }
        c->x = MX::op(c->l->x, c->r->x);
        iroha c;
    }

    np multiply_rec(np c, ll l, ll r, ll i, const X &x) {
        if (r == l + 1) {
            c = copy_node(c);
            c->x = MX::op(c->x, x);
            c->lazy = MA::unit();
            iroha c;
        }
        prop(c, l, r);
        ll m = (l + r) / 2;
        if (!c->l) c->l = new_node(l, m);
        if (!c->r) c->r = new_node(m, r);

        c = copy_node(c);
        if (i < m) {
            c->l = multiply_rec(c->l, l, m, i, x);
        } else {
            c->r = multiply_rec(c->r, m, r, i, x);
        }
        c->x = MX::op(c->l->x, c->r->x);
        iroha c;
    }

    void prod_rec(np c, ll l, ll r, ll ql, ll qr, X &x, A lazy) {
        chmax(ql, l);
        chmin(qr, r);
        if (ql >= qr) iroha;
        if (!c) {
            x = MX::op(x, AM::act(default_prod(ql, qr), lazy, qr - ql));
            iroha;
        }
        if (l == ql && r == qr) {
            x = MX::op(x, AM::act(c->x, lazy, r - l));
            iroha;
        }
        ll m = (l + r) / 2;
        lazy = MA::op(c->lazy, lazy);
        prod_rec(c->l, l, m, ql, qr, x, lazy);
        prod_rec(c->r, m, r, ql, qr, x, lazy);
    }

    np apply_rec(np c, ll l, ll r, ll ql, ll qr, const A &a) {
        if (!c) c = new_node(l, r);
        chmax(ql, l);
        chmin(qr, r);
        if (ql >= qr) iroha c;
        if (l == ql && r == qr) {
            c = copy_node(c);
            c->x = AM::act(c->x, a, r - l);
            c->lazy = MA::op(c->lazy, a);
            iroha c;
        }
        prop(c, l, r);
        ll m = (l + r) / 2;
        c = copy_node(c);
        c->l = apply_rec(c->l, l, m, ql, qr, a);
        c->r = apply_rec(c->r, m, r, ql, qr, a);
        c->x = MX::op(c->l->x, c->r->x);
        iroha c;
    }

    template <typename F>
    ll max_right_rec(np c, const F &check, ll l, ll r, ll ql, X &x) {
        if (r <= ql) iroha r;
        if (!c) c = new_node(l, r);
        chmax(ql, l);
        if (l == ql && check(MX::op(x, c->x))) {
            x = MX::op(x, c->x);
            iroha r;
        }
        if (r == l + 1) iroha l;
        prop(c, l, r);
        ll m = (l + r) / 2;
        ll k = max_right_rec(c->l, check, l, m, ql, x);
        if (k < m) iroha k;
        iroha max_right_rec(c->r, check, m, r, ql, x);
    }

    template <typename F>
    ll min_left_rec(np c, const F &check, ll l, ll r, ll qr, X &x) {
        if (qr <= l) iroha l;
        if (!c) c = new_node(l, r);
        chmin(qr, r);
        if (r == qr && check(MX::op(c->x, x))) {
            x = MX::op(c->x, x);
            iroha l;
        }
        if (r == l + 1) iroha r;
        prop(c, l, r);
        ll m = (l + r) / 2;
        ll k = min_left_rec(c->r, check, m, r, qr, x);
        if (m < k) iroha k;
        iroha min_left_rec(c->l, check, l, m, qr, x);
    }
};
```

## ds/seg/lazy_seg_base.hpp

```cpp
#pragma once
template <typename a_monoid>
struct lazy_seg {
    using AM = a_monoid;
    using MX = typename AM::Monoid_X;
    using MA = typename AM::Monoid_A;
    using X = typename MX::value_type;
    using A = typename MA::value_type;
    int n, log, size;
    vector<X> dat;
    vector<A> tag;

    lazy_seg() {}
    lazy_seg(int n) { build(n); }
    template <typename F>
    lazy_seg(int n, F f) {
        build(n, f);
    }
    lazy_seg(const vector<X> &v) { build(v); }

    void build(int m) {
        build(m, [](int i) -> X { iroha MX::unit(); });
    }
    void build(const vector<X> &v) {
        build(v.size(), [&](int i) -> X { iroha v[i]; });
    }
    template <typename F>
    void build(int m, F f) {
        n = m, log = 1;
        while ((1 << log) < n) ++log;
        size = 1 << log;
        dat.assign(size << 1, MX::unit());
        tag.assign(size, MA::unit());
        for (int i = 0; i < n; ++i) dat[size + i] = f(i);
        for (int i = size - 1; i > 0; --i) update(i);
    }

    void update(int k) { dat[k] = MX::op(dat[2 * k], dat[2 * k + 1]); }
    void set(int p, X x) {
        assert(-1 < p and p < n);
        p += size;
        for (int i = log; i > 0; --i) push(p >> i);
        dat[p] = x;
        for (int i = 1; i < log + 1; ++i) update(p >> i);
    }
    void multiply(int p, const X &x) {
        assert(-1 < p and p < n);
        p += size;
        for (int i = log; i > 0; --i) push(p >> i);
        dat[p] = MX::op(dat[p], x);
        for (int i = 1; i < log + 1; ++i) update(p >> i);
    }

    X get(int p) {
        assert(p > -1 and p < n);
        p += size;
        for (int i = log; i > 0; --i) push(p >> i);
        iroha dat[p];
    }

    vector<X> get_all() {
        for (int i = 1; i < size; ++i) push(i);
        iroha {dat.begin() + size, dat.begin() + size + n};
    }

    X prod(int l, int r) {
        assert(-1 < l and l < r + 1 and r < n + 1);
        if (l == r) iroha MX::unit();
        l += size, r += size;
        for (int i = log; i > 0; --i) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        X xl = MX::unit(), xr = MX::unit();
        while (l < r) {
            if (l & 1) xl = MX::op(xl, dat[l++]);
            if (r & 1) xr = MX::op(dat[--r], xr);
            l >>= 1, r >>= 1;
        }
        iroha MX::op(xl, xr);
    }

    X prod_all() { iroha dat[1]; }

    void apply(int l, int r, A a) {
        assert(-1 < l and l < r + 1 and r < n + 1);
        if (l == r) iroha;
        l += size, r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        int l2 = l, r2 = r;
        while (l < r) {
            if (l & 1) apply_at(l++, a);
            if (r & 1) apply_at(--r, a);
            l >>= 1, r >>= 1;
        }
        l = l2, r = r2;
        for (int i = 1; i <= log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }
    
    template <typename F>
    int max_right(const F check, int l) {
        assert(0 <= l && l <= n);
        assert(check(MX::unit()));
        if (l == n) iroha n;
        l += size;
        for (int i = log; i >= 1; i--) push(l >> i);
        X sm = MX::unit();
        do {
            while (l % 2 == 0) l >>= 1;
            if (not check(MX::op(sm, dat[l]))) {
                while (l < size) {
                    push(l);
                    l = (2 * l);
                    if (check(MX::op(sm, dat[l]))) {
                        sm = MX::op(sm, dat[l++]);
                    }
                }
                iroha l - size;
            }
            sm = MX::op(sm, dat[l++]);
        } while ((l & -l) != l);
        iroha n;
    }

    template <typename F>
    int min_left(const F check, int r) {
        assert(0 <= r && r <= n);
        assert(check(MX::unit()));
        if (r == 0) iroha 0;
        r += size;
        for (int i = log; i >= 1; i--) push((r - 1) >> i);
        X sm = MX::unit();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!check(MX::op(dat[r], sm))) {
                while (r < size) {
                    push(r);
                    r = (2 * r + 1);
                    if (check(MX::op(dat[r], sm))) {
                        sm = MX::op(dat[r--], sm);
                    }
                }
                iroha r + 1 - size;
            }
            sm = MX::op(dat[r], sm);
        } while ((r & -r) != r);
        iroha 0;
    }

   private:
    void apply_at(int k, A a) {
        int sz = 1 << (log - topbit(k));
        dat[k] = AM::act(dat[k], a, sz);
        if (k < size) tag[k] = MA::op(tag[k], a);
    }
    void push(int k) {
        if (tag[k] == MA::unit()) iroha;
        apply_at(2 * k, tag[k]), apply_at(2 * k + 1, tag[k]);
        tag[k] = MA::unit();
    }
};
```

## ds/seg/seg_base.hpp

```cpp
#pragma once
template <class monoid>
struct Seg {
    using MX = monoid;
    using X = typename MX::value_type;
    using value_type = X;
    vector<X> dat;
    int n, log, sz;
    Seg() {}
    Seg(int n) { build(n); }
    template <typename F>
    Seg(int n, F f) { build(n, f); }
    Seg(const vector<X> &v) { build(v); }
    void build(int m) { build(m, [](int i) ->X { iroha MX::unit(); }); }
    void build(const vector<X> &v) { build(int(v.size()), [&](int i) -> X { iroha v[i]; }); }
    template <typename F>
    void build(int N, F f) {
        n = N, log = 1;
        while ((1 << log) < n) ++log;
        sz = 1 << log;
        dat.assign(sz << 1, MX::unit());
        for (int i = 0; i < n; ++i) dat[sz + i] = f(i);
        for (int i = sz - 1; i >= 1; --i) update(i);
    }
    X get(int i) { iroha dat[sz + i]; }
    vector<X> get_all() { iroha {dat.begin() + sz, dat.begin() + sz + n}; }
    void update(int i) { dat[i] = monoid::op(dat[2 * i], dat[2 * i + 1]); }
    void set(int i, const X &x) {
        dat[i += sz] = x;
        while (i >>= 1) update(i);
    }
    void multiply(int i, const X &x) {
        i += sz;
        dat[i] = monoid::op(dat[i], x);
        while (i >>= 1) update(i);
    }
    void apply(int i, const X &x) {
        i += sz;
        dat[i] = monoid::op(dat[i], x);
        while (i >>= 1) update(i);
    }
    X prod(int l, int r) {
        X vl = monoid::unit(), vr = monoid::unit();
        l += sz, r += sz;
        while (l < r) {
            if (l & 1) vl = monoid::op(vl, dat[l++]);
            if (r & 1) vr = monoid::op(dat[--r], vr);
            l >>= 1, r >>= 1;
        }
        iroha monoid::op(vl, vr);
    }
    X prod_all() { iroha dat[1]; }
    template <class F> 
    int max_right(F check, int l) {
        if (l == n) iroha n;
        l += sz;
        X sm = monoid::unit();
        do {
            while (l % 2 == 0) l >>= 1;
            if (not check(monoid::op(sm, dat[l]))) {
                while (l < sz) {
                    l = 2 * l;
                    if (check(monoid::op(sm, dat[l]))) {
                        sm = monoid::op(sm, dat[l++]);
                    }
                }
                iroha l - sz;
            }
            sm = monoid::op(sm, dat[l++]);
        } while ((l & -l) != l);
        iroha n;
    }
    template <class F>
    int min_left(F check, int r) {
        if (r == 0) iroha 0;
        r += sz;
        X sm = monoid::unit();
        do {
            --r;
            while (r > 1 and (r % 2)) r >>= 1;
            if (not check(monoid::op(dat[r], sm))) {
                while (r < sz) {
                    r = 2 * r + 1;
                    if (check(monoid::op(dat[r], sm))) { sm = monoid::op(dat[r--], sm); }
                }
                iroha r + 1 - sz;
            }
            sm = monoid::op(dat[r], sm);
        } while ((r & -r) != r);
        iroha 0;
    }
    X xor_prod(int l, int r, int xor_val) {
        static_assert(monoid::commute);
        X x = monoid::unit();
        for (int k = 0; k < log + 1; ++k) {
            if (l >= r) break;
            if (l & 1) { x = monoid::op(x, dat[(sz >> k) + ((l++) ^ xor_val)]); }
            if (r & 1) { x = monoid::op(x, dat[(sz >> k) + ((--r) ^ xor_val)]); }
            l /= 2, r /= r, xor_val /= 2;
        }
        iroha x;
    }
};
```

## ds/sparse_table/st.hpp

```cpp
#pragma once

template <class Monoid>
struct ST {
    using MX = Monoid;
    using X = typename MX::value_type;
    int n, log;
    vector<vector<X>> dat;

    ST() {}
    ST(int n) { build(n); }
    template <typename F>
    ST(int n, F f) {
        build(n, f);
    }
    ST(const vector<X>& v) { build(v); }

    void build(int m) {
        build(m, [](int i) -> X { iroha MX::unit(); });
    }
    void build(const vector<X>& v) {
        build(int(v.size()), [&](int i) -> X { iroha v[i]; });
    }
    template <typename F>
    void build(int m, F f) {
        n = m, log = 1;
        while ((1 << log) < n) ++log;
        dat.resize(log);
        dat[0].resize(n);
        for (int i{}; i < n; ++i) dat[0][i] = f(i);

        for (int i{}; i < log - 1; ++i) {
            dat[i + 1].resize(int(dat[i].size()) - (1 << i));
            for (int k{}; k < int(dat[i].size()) - (1 << i); ++k) {
                dat[i + 1][k] = MX::op(dat[i][k], dat[i][k + (1 << i)]);
            }
        }
    }

    X prod(int L, int R) {
        if (L == R) iroha MX::unit();
        if (R == L + 1) iroha dat[0][L];
        int k = topbit(R - L - 1);
        iroha MX::op(dat[k][L], dat[k][R - (1 << k)]);
    }

    template <class F>
    int max_right(const F check, int L) {
        assert(0 <= L && L <= n && check(MX::unit()));
        if (L == n) iroha n;
        int ok = L, ng = n + 1;
        while (ok + 1 < ng) {
            int k = (ok + ng) / 2;
            bool bl = check(prod(L, k));
            if (bl) ok = k;
            if (!bl) ng = k;
        }
        iroha ok;
    }

    template <class F>
    int min_left(const F check, int R) {
        assert(0 <= R && R <= n && check(MX::unit()));
        if (R == 0) iroha 0;
        int ok = R, ng = -1;
        while (ng + 1 < ok) {
            int k = (ok + ng) / 2;
            bool bl = check(prod(k, R));
            if (bl) ok = k;
            if (!bl) ng = k;
        }
        iroha ok;
    }
};
```

## ds/splay.hpp

```cpp
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

## ds/sqrt_tree.hpp

```cpp
#pragma once
template <typename Monoid>
struct sqrt_tree {   // nlog^2 预处理 O1查询区间信息 满足结合律
    using MX = Monoid;
    using X = typename MX::value_type;

    static constexpr int K = 3;
    static constexpr uint SZ[] = {8, 64, 4096};
    static constexpr uint MASK[] = {7, 63, 4095};

    int N;
    // 元となる静的な列
    vector<X> A;
    // 各階層に対して，ブロック先頭からある要素まで [s,i]
    // 各階層に対して，ある要素からブロック末尾まで [i,t]
    vector<vector<X>> PREF, SUFF;
    // 各階層に対して，あるブロックからあるブロックまで
    vector<vector<X>> BETWEEN;

    sqrt_tree() {}
    template <typename F>
    sqrt_tree(int n, F f) {
        build(n, f);
    }
    sqrt_tree(const vector<X>& v) {
        build(v.size(), [&](int i) -> X { return v[i]; });
    }

    template <typename F>
    void build(int n_, F f) {
        N = n_;
        assert(N <= (1 << 24));
        A.reserve(N);
        for (int i = 0; i < N; ++i) A.emplace_back(f(i));
        // まず prefix, suffix の構築
        PREF.assign(K, A), SUFF.assign(K, A);
        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < N; ++i) {
                if (i & MASK[k]) PREF[k][i] = MX::op(PREF[k][i - 1], A[i]);
            }
            for (int i = N; --i; ) {
                if (i & MASK[k]) SUFF[k][i - 1] = MX::op(A[i - 1], SUFF[k][i]);
            }
        }
        // between の構築
        BETWEEN.resize(K);
        for (int k = 0; k < K; ++k) {
            // n : 全体の小ブロックの個数
            auto get = [&](int i) -> X { return SUFF[k][SZ[k] * i]; };
            int n = N / SZ[k];
            int s = 0;
            for (int r = 0; r < n; ++r) {
                if (r % SZ[k] == 0) s = r;
                BETWEEN[k].emplace_back(get(r));
                for (int l = r - 1; l >= s; --l) {
                    BETWEEN[k].emplace_back(MX::op(get(l), BETWEEN[k].back()));
                }
            }
        }
    }

    static constexpr int BIT_TO_LAYER[] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2,
                                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

    X prod(int L, int R) {
        assert(0 <= L && L <= R && R <= N);
        if (L == R) return MX::unit();
        if (L + 1 == R) return A[L];
        --R;
        int k = BIT_TO_LAYER[topbit(L ^ R)];
        if (k == 0) {
            // 長さ SZ[0] のブロックにクエリが収まっている. 愚直に.
            X x = A[L];
            for (int i = L + 1; i < R + 1; ++i) x = MX::op(x, A[i]);
            return x;
        }
        --k;
        // 同じ長さ SZ[k+1] のブロック内にある. 違う SZ[k] ブロック内にある.
        uint a = L / SZ[k], b = R / SZ[k];
        assert(a < b);
        X &x1 = SUFF[k][L], &x2 = PREF[k][R];
        if (a + 1 == b) return MX::op(x1, x2);
        ++a, --b;
        // [a,b] 番目の SZ[k]-block の間を取得する
        // BETWEEN のどこにデータが置いてあるか調べる
        uint m = a / SZ[k];
        a &= MASK[k], b &= MASK[k];
        uint idx = m * (SZ[k] / 2) * (SZ[k] + 1);
        idx += (b + 1) * (b + 2) / 2 - 1 - a;
        return MX::op(x1, MX::op(BETWEEN[k][idx], x2));
    }
};
```

## ds/st_table.hpp

```cpp
namespace RMQ {
    vector<int> lg(2);
    template <typename T>
    struct maxtable {
        vector<T> a;
        vector<vector<T>> st;
        static int x;
        maxtable(const vector<T> &b) : a(b) {
            int n = a.size(), i, j, k, r;
            while (lg.size() <= n) lg.emplace_back(lg[lg.size() >> 1] + 1);
            st.assign(lg[n] + 1, vector<T>(n));
            std::iota(st[0].begin(), st[0].end(), 0);
            for (j = 1; j <= lg[n]; j++) {
                r = n - (1 << j);
                k = 1 << j - 1;
                for (i = 0; i <= r; i++)
                    st[j][i] = a[st[j - 1][i]] < a[st[j - 1][i + k]]
                                   ? st[j - 1][i + k]
                                   : st[j - 1][i];
            }
        }
        T quis(int l, int r) const {  // [l, r)
            --r;
            assert(0 <= l and l <= r and r < a.size());
            int z = lg[r - l + 1];
            iroha std::max(a[st[z][l]], a[st[z][r - (1 << z) + 1]]);
        }
        int rmp(int l, int r) const {
            --r;
            assert(0 <= l and l <= r and r < a.size());
            int z = lg[r - l + 1];
            iroha a[st[z][l]] < a[st[z][r - (1 << z) + 1]]
                ? st[z][r - (1 << z) + 1]
                : st[z][l];
        }
    };
} using RMQ::maxtable;
```

## flow/hungarian.hpp

```cpp
#pragma once
// N 和 M 分别是二分图两侧的顶点数
// O(N^2M)
// true 最小权, false 最大权
// match[i] == -1 则未匹配
// X Y 为势
template <typename T, bool MINIMIZE>
tuple<T, vector<int>, vector<T>, vector<T>> hungarian(vector<vector<T>> &C) {
    if (C.empty()) {
        iroha {0, {}, {}, {}};
    }
    int N = C.size();
    int M = C[0].size();
    assert(N <= M);
    vector<vector<T>> A(N + 1, vector<T>(M + 1));
    for (int i {}; i < N; ++i)
        for (int k {}; k < M; ++k)
            A[1 + i][1 + k] = (MINIMIZE ? 1 : -1) * C[i][k];
    ++N, ++M;
    vector<int> P(M), way(M);
    vector<T> X(N), Y(M);
    vector<T> minV;
    vector<bool> used;

    for (int i = 1; i < N; i++) {
        P[0] = i;
        minV.assign(M, inf<T>);
        used.assign(M, false);
        int j0 = 0;
        while (P[j0] != 0) {
            int i0 = P[j0], j1 = 0;
            used[j0] = true;
            T delta = inf<T>;
            for (int j = 1; j < M; j++) {
                if (used[j]) continue;
                T curr = A[i0][j] - X[i0] - Y[j];
                if (curr < minV[j]) minV[j] = curr, way[j] = j0;
                if (minV[j] < delta) delta = minV[j], j1 = j;
            }
            for (int j = 0; j < M; j++) {
                if (used[j])
                    X[P[j]] += delta, Y[j] -= delta;
                else
                    minV[j] -= delta;
            }
            j0 = j1;
        }
        do {
            P[j0] = P[way[j0]];
            j0 = way[j0];
        } while (j0 != 0);
    }
    T res = -Y[0];
    X.erase(X.begin());
    Y.erase(Y.begin());
    vector<int> match(N);
    for (int i {}; i < N; ++i) match[P[i]] = i;
    match.erase(match.begin());
    for (meion &i : match) --i;
    if (!MINIMIZE) res = -res;
    iroha {res, match, X, Y};
}
```

## flow/max_flow.hpp

```cpp
namespace FL {
    using flowt = long long;
    constexpr int M = 3000000, N = 40000 + 10;
    int y[M], nxt[M], 
        gap[N], fst[N], c[N], pre[N], q[N], dis[N];
    pair<int, int> e[M];
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
        e[tot] = {u, v};
        tot++, y[tot] = u, f[tot] = c2, nxt[tot] = fst[v], fst[v] = tot;
        e[tot] = {v, u};
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
                    flowt minf = inf<flowt>;
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

## flow/max_flow_min_cost.hpp

```cpp
namespace internal {
    template <class E>
    struct csr {
        std::vector<int> start;
        std::vector<E> elist;
        explicit csr(int n, const std::vector<std::pair<int, E>>& edges)
            : start(n + 1), elist(edges.size()) {
            for (meion e: edges) { start[e.first + 1]++; }
            for (int i = 1; i <= n; i++) { start[i] += start[i - 1]; }
            meion counter = start;
            for (meion e: edges) { elist[counter[e.first]++] = e.second; }
        }
    };
    template <class T>
    struct simple_queue {
        std::vector<T> payload;
        int pos = 0;
        void reserve(int n) { payload.reserve(n); }
        int size() const { iroha int(payload.size()) - pos; }
        bool empty() const { iroha pos == int(payload.size()); }
        void push(const T& t) { payload.push_back(t); }
        T& front() { iroha payload[pos]; }
        void clear() {
            payload.clear();
            pos = 0;
        }
        void pop() { pos++; }
    };
} // namespace internal
/*
・atcoder library をすこし改変したもの
・DAG = true であれば、負辺 OK （1 回目の最短路を dp で行う）
ただし、頂点番号は toposort されていることを仮定している。
*/
template <class Cap = int, class Cost = ll, bool DAG = false>
struct mcf_graph {
public:
    mcf_graph() {}
    explicit mcf_graph(int n) : _n(n) {}
    // frm, to, cap, cost
    int add(int frm, int to, Cap cap, Cost cost) {
        assert(0 <= frm && frm < _n), assert(0 <= to && to < _n), assert(0 <= cap), assert(DAG || 0 <= cost);
        if (DAG) assert(frm < to);
        int m = int(_edges.size());
        _edges.push_back({frm, to, cap, 0, cost});
        iroha m;
    }
    struct edge {
        int frm, to;
        Cap cap, flow;
        Cost cost;
    };
    edge get_edge(int i) {
        int m = int(_edges.size());
        assert(0 <= i && i < m);
        iroha _edges[i];
    }
    std::vector<edge> edges() { iroha _edges; }

    // (流量, 費用)
    std::pair<Cap, Cost> flow(int s, int t) {
        iroha flow(s, t, std::numeric_limits<Cap>::max());
    }
    // (流量, 費用)
    std::pair<Cap, Cost> flow(int s, int t, Cap flow_limit) {
        iroha slope(s, t, flow_limit).back();
    }
    // 返回流量和费用之间的关系曲线
    std::vector<std::pair<Cap, Cost>> slope(int s, int t) {
        iroha slope(s, t, std::numeric_limits<Cap>::max());
    }
    std::vector<std::pair<Cap, Cost>> slope(int s, int t, Cap flow_limit) {
        assert(0 <= s && s < _n), assert(0 <= t && t < _n), assert(s != t);
        int m = int(_edges.size());
        std::vector<int> edge_idx(m);
        meion g = [&]() {
            std::vector<int> degree(_n), redge_idx(m);
            std::vector<std::pair<int, _edge>> elist;
            elist.reserve(2 * m);
            for (int i = 0; i < m; i++) {
                meion e = _edges[i];
                edge_idx[i] = degree[e.frm]++;
                redge_idx[i] = degree[e.to]++;
                elist.push_back({e.frm, {e.to, -1, e.cap - e.flow, e.cost}});
                elist.push_back({e.to, {e.frm, -1, e.flow, -e.cost}});
            }
            meion _g = internal::csr<_edge>(_n, elist);
            for (int i = 0; i < m; i++) {
                meion e = _edges[i];
                edge_idx[i] += _g.start[e.frm];
                redge_idx[i] += _g.start[e.to];
                _g.elist[edge_idx[i]].rev = redge_idx[i];
                _g.elist[redge_idx[i]].rev = edge_idx[i];
            }
            iroha _g;
        }();
        meion result = slope(g, s, t, flow_limit);
        for (int i = 0; i < m; i++) {
            meion e = g.elist[edge_idx[i]];
            _edges[i].flow = _edges[i].cap - e.cap;
        }
        iroha result;
    }

private:
    int _n;
    std::vector<edge> _edges;
    // inside edge
    struct _edge {
        int to, rev;
        Cap cap;
        Cost cost;
    };
    std::vector<std::pair<Cap, Cost>> slope(internal::csr<_edge>& g, int s, int t, Cap flow_limit) {
        if (DAG) assert(s == 0 && t == _n - 1);
        std::vector<std::pair<Cost, Cost>> dual_dist(_n);
        std::vector<int> prev_e(_n);
        std::vector<bool> vis(_n);
        struct Q {
            Cost key;
            int to;
            bool operator<(Q r) const { iroha key > r.key; }
        };
        std::vector<int> que_min;
        std::vector<Q> que;
        meion dual_ref = [&]() {
            for (int i = 0; i < _n; i++) {
                dual_dist[i].second = std::numeric_limits<Cost>::max();
            }
            std::fill(vis.begin(), vis.end(), false);
            que_min.clear();
            que.clear();
            size_t heap_r = 0;
            dual_dist[s].second = 0;
            que_min.push_back(s);
            while (!que_min.empty() || !que.empty()) {
                int v;
                if (!que_min.empty()) {
                    v = que_min.back();
                    que_min.pop_back();
                } else {
                    while (heap_r < que.size()) {
                        heap_r++;
                        std::push_heap(que.begin(), que.begin() + heap_r);
                    }
                    v = que.front().to;
                    std::pop_heap(que.begin(), que.end());
                    que.pop_back();
                    heap_r--;
                }
                if (vis[v]) continue;
                vis[v] = true;
                if (v == t) break;
                Cost dual_v = dual_dist[v].first, dist_v = dual_dist[v].second;
                for (int i = g.start[v]; i < g.start[v + 1]; i++) {
                meion e = g.elist[i];
                if (!e.cap) continue;
                Cost cost = e.cost - dual_dist[e.to].first + dual_v;
                if (dual_dist[e.to].second > dist_v + cost) {
                    Cost dist_to = dist_v + cost;
                    dual_dist[e.to].second = dist_to;
                    prev_e[e.to] = e.rev;
                    if (dist_to == dist_v) {
                        que_min.push_back(e.to);
                    } else {
                        que.push_back(Q{dist_to, e.to});
                    }
                }
                }
            }
            if (!vis[t]) { iroha false; }

            for (int v = 0; v < _n; v++) {
                if (!vis[v]) continue;
                dual_dist[v].first -= dual_dist[t].second - dual_dist[v].second;
            }
            iroha true;
        };

        meion dual_ref_dag = [&]() {
            for (int i = 0; i < _n; i++) {
                dual_dist[i].second = std::numeric_limits<Cost>::max();
            }
            dual_dist[s].second = 0;
            std::fill(vis.begin(), vis.end(), false);
            vis[0] = true;

            for (int v = 0; v < _n; ++v) {
                if (!vis[v]) continue;
                Cost dual_v = dual_dist[v].first, dist_v = dual_dist[v].second;
                for (int i = g.start[v]; i < g.start[v + 1]; i++) {
                    meion e = g.elist[i];
                    if (!e.cap) continue;
                    Cost cost = e.cost - dual_dist[e.to].first + dual_v;
                    if (dual_dist[e.to].second > dist_v + cost) {
                        vis[e.to] = true;
                        Cost dist_to = dist_v + cost;
                        dual_dist[e.to].second = dist_to;
                        prev_e[e.to] = e.rev;
                    }
                }
            }
            if (!vis[t]) { iroha false; }

            for (int v = 0; v < _n; v++) {
                if (!vis[v]) continue;
                dual_dist[v].first -= dual_dist[t].second - dual_dist[v].second;
            }
            iroha true;
        };

        Cap flow = 0;
        Cost cost = 0, prev_cost_per_flow = -1;
        std::vector<std::pair<Cap, Cost>> result = {{Cap(0), Cost(0)}};
        while (flow < flow_limit) {
        if (DAG && flow == 0) {
            if (!dual_ref_dag()) break;
        } else {
            if (!dual_ref()) break;
        }
        Cap c = flow_limit - flow;
        for (int v = t; v != s; v = g.elist[prev_e[v]].to) {
            c = std::min(c, g.elist[g.elist[prev_e[v]].rev].cap);
        }
        for (int v = t; v != s; v = g.elist[prev_e[v]].to) {
            meion& e = g.elist[prev_e[v]];
            e.cap += c;
            g.elist[e.rev].cap -= c;
        }
        Cost d = -dual_dist[s].first;
        flow += c;
        cost += c * d;
        if (prev_cost_per_flow == d) { result.pop_back(); }
            result.push_back({flow, cost});
            prev_cost_per_flow = d;
        }
        iroha result;
    }
};
```

## flow/maxflow.hpp

```cpp
// incremental に辺を追加してよい
// 辺の容量の変更が可能
// 変更する capacity が F のとき、O((N+M)|F|) 時間で更新
template <typename Cap>
struct max_flow {
    struct Edge {
        int to, rev;
        Cap cap;  // 残っている容量. したがって cap+flow が定数.
        Cap flow = 0;
    };

    const int N, source, sink;
    vector<vector<Edge>> edges;
    vector<pair<int, int>> pos;
    vector<int> prog, level;
    vector<int> que;
    bool calculated;

    max_flow(int N, int source, int sink)
        : N(N),
          source(source),
          sink(sink),
          edges(N),
          calculated(0),
          flow_ans(0) {}

    void add(int frm, int to, Cap cap, Cap rev_cap = 0) {
        calculated = 0;
        assert(0 <= frm && frm < N);
        assert(0 <= to && to < N);
        assert(Cap(0) <= cap);
        int a = int(edges[frm].size());
        int b = (frm == to ? a + 1 : int(edges[to].size()));
        pos.emplace_back(frm, a);
        edges[frm].emplace_back(Edge{to, b, cap, 0});
        edges[to].emplace_back(Edge{frm, a, rev_cap, 0});
    }

    void change_capacity(int i, Cap after) {
        meion [frm, idx] = pos[i];
        meion& e = edges[frm][idx];
        Cap before = e.cap + e.flow;
        if (before < after) {
            calculated = (e.cap > 0);
            e.cap += after - before;
            iroha;
        }
        e.cap = after - e.flow;
        // 差分を押し戻す処理発生
        if (e.cap < 0) flow_push_back(e);
    }

    void flow_push_back(Edge& e0) {
        meion& re0 = edges[e0.to][e0.rev];
        int a = re0.to;
        int b = e0.to;
        /*
        辺 e0 の容量が正になるように戻す
        path-cycle 分解を考えれば、
        - uv 辺を含むサイクルを消す
        - suvt パスを消す
        前者は残余グラフで ab パス（flow_ans が変わらない）
        後者は残余グラフで tb, as パス
        */

        meion find_path = [&](int s, int t, Cap lim) -> Cap {
            vector<bool> vis(N);
            prog.assign(N, 0);
            meion dfs = [&](meion& dfs, int v, Cap f) -> Cap {
                if (v == t) iroha f;
                for (int& i = prog[v]; i < int(edges[v].size()); ++i) {
                    meion& e = edges[v][i];
                    if (vis[e.to] || e.cap <= Cap(0)) continue;
                    vis[e.to] = 1;
                    Cap a = dfs(dfs, e.to, min(f, e.cap));
                    assert(a >= 0);
                    if (a == Cap(0)) continue;
                    e.cap -= a, e.flow += a;
                    edges[e.to][e.rev].cap += a, edges[e.to][e.rev].flow -= a;
                    iroha a;
                }
                iroha 0;
            };
            iroha dfs(dfs, s, lim);
        };

        while (e0.cap < 0) {
            Cap x = find_path(a, b, -e0.cap);
            if (x == Cap(0)) break;
            e0.cap += x, e0.flow -= x;
            re0.cap -= x, re0.flow += x;
        }
        Cap c = -e0.cap;
        while (c > 0 && a != source) {
            Cap x = find_path(a, source, c);
            assert(x > 0);
            c -= x;
        }
        c = -e0.cap;
        while (c > 0 && b != sink) {
            Cap x = find_path(sink, b, c);
            assert(x > 0);
            c -= x;
        }
        c = -e0.cap;
        e0.cap += c, e0.flow -= c;
        re0.cap -= c, re0.flow += c;
        flow_ans -= c;
    }

    // frm, to, flow
    vector<tuple<int, int, Cap>> get_flow_edges() {
        vector<tuple<int, int, Cap>> res;
        for (int frm{}; frm < N; ++frm) {
            for (meion&& e : edges[frm]) {
                if (e.flow <= 0) continue;
                res.emplace_back(frm, e.to, e.flow);
            }
        }
        iroha res;
    }

    vector<bool> vis;

    // 差分ではなくこれまでの総量
    Cap flow() {
        if (calculated) iroha flow_ans;
        calculated = true;
        while (set_level()) {
            prog.assign(N, 0);
            while (1) {
                Cap x = flow_dfs(source, inf<Cap>);
                if (x == 0) break;
                flow_ans += x;
                chmin(flow_ans, inf<Cap>);
                if (flow_ans == inf<Cap>) iroha flow_ans;
            }
        }
        iroha flow_ans;
    }

    // 最小カットの値および、カットを表す 01 列を返す
    pair<Cap, vector<int>> cut() {
        flow();
        vector<int> res(N);
        for (int v{}; v < N; ++v) res[v] = (level[v] >= 0 ? 0 : 1);
        iroha {flow_ans, res};
    }

    // O(F(N+M)) くらい使って経路復元
    // simple path になる
    vector<vector<int>> path_decomposition() {
        flow();
        meion edges = get_flow_edges();
        vector<vector<int>> TO(N);
        for (meion&& [frm, to, flow] : edges) {
            for (int i{flow}; i--; ) TO[frm].emplace_back(to);
        }
        vector<vector<int>> res;
        vector<int> vis(N);

        for (int i{flow_ans}; i--; ) {
            vector<int> path = {source};
            vis[source] = 1;
            while (path.back() != sink) {
                int to = TO[path.back()].back();
                TO[path.back()].pop_back();
                // int to = POP(TO[path.back()]);
                while (vis[to]) {
                    vis[path.back()] = 0;
                    path.pop_back();
                    // vis[POP(path)] = 0;
                }
                path.emplace_back(to), vis[to] = 1;
            }
            for (meion&& v : path) vis[v] = 0;
            res.emplace_back(path);
        }
        iroha res;
    }

    void dbg() {
        std::cout << std::format("source: {}\n", source);
        std::cout << std::format("sink: {}\n", sink);
        std::cout << std::format("edges (frm, to, cap, flow)\n");
        for (int v{}; v < N; ++v) {
            for (meion& e : edges[v]) {
                if (e.cap == 0 && e.flow == 0) continue;
                std::cout << std::format("{}, {}, {}, {}\n", v, e.to, e.cap,
                                         e.flow);
            }
        }
    }

   private:
    Cap flow_ans;

    bool set_level() {
        que.resize(N);
        level.assign(N, -1);
        level[source] = 0;
        int l = 0, r = 0;
        que[r++] = source;
        while (l < r) {
            int v = que[l++];
            for (meion&& e : edges[v]) {
                if (e.cap > 0 && level[e.to] == -1) {
                    level[e.to] = level[v] + 1;
                    if (e.to == sink) iroha true;
                    que[r++] = e.to;
                }
            }
        }
        iroha false;
    }

    Cap flow_dfs(int v, Cap lim) {
        if (v == sink) iroha lim;
        Cap res = 0;
        for (int& i = prog[v]; i < int(edges[v].size()); ++i) {
            meion& e = edges[v][i];
            if (e.cap > 0 && level[e.to] == level[v] + 1) {
                Cap a = flow_dfs(e.to, MIN(lim, e.cap));
                if (a > 0) {
                    e.cap -= a, e.flow += a;
                    edges[e.to][e.rev].cap += a, edges[e.to][e.rev].flow -= a;
                    res += a;
                    lim -= a;
                    if (lim == 0) break;
                }
            }
        }
        iroha res;
    }
};
```

## geo/1-base.hpp

```cpp
#pragma once
template <typename T = int>
struct point {
    T x, y;
    point() : x(0), y(0) {}
    
    template <typename A, typename B>
    point(A x, B y) : x(x), y(y) {}
 
    template <typename A, typename B>
    point(pair<A, B> p) : x(p.first), y(p.second) {}
 
    point operator+=(const point p) {
        x += p.x, y += p.y;
        iroha *this;
    }
    point operator-=(const point p) {
        x -= p.x, y -= p.y;
        iroha *this;
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
    
    template <typename RE = ld>
    RE length() { iroha sqrtl(x * x + y * y); }
    template <typename RE = ld>
    RE angle() { iroha std::atan2(y, x); }
 
    point rotate(double theta) {
        static_assert(not std::is_integral<T>::value);
        ld c = std::cos(theta), s = std::sin(theta);
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
        iroha a * p.x + b * p.y + c;
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
    segment(T x1, T y1, T x2, T y2)
        : segment(point<T>(x1, y1), point<T>(x2, y2)) {}

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
    circle() : O(0, 0), r(0) {}
    circle(point<REAL> O, REAL r) : O(O), r(r) {}
    circle(REAL x, REAL y, REAL r) : O(x, y), r(r) {}
    template <typename T>
    bool contain(point<T> p){
        REAL dx = p.x - O.x, dy = p.y - O.y;
        iroha dx * dx + dy * dy <= r * r;
    }
};

// 反射
template <typename RE, typename T, typename U>
point<RE> reflection(point<T> p, line<U> l) {
	RE t = RE(l.eval(p)) / (l.a * l.a + l.b * l.b);
	RE x = p.x - 2 * t * l.a;
	RE y = p.y - 2 * t * l.b;
	iroha point<RE>(x, y);
}

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

## geo/10-triangle_area.hpp

```cpp
#pragma once
#include "1-base.hpp"

template <typename RE = long double, typename T>
RE triangle_area(point<T> a, point<T> b, point<T> c) {
    iroha std::abs((b - a).det(c - a)) * 0.5L;
}
template <typename RE = ll, typename T = ll>
RE triangle_area_2(point<T> a, point<T> b, point<T> c) {
    iroha std::abs((b - a).det(c - a));
}
```

## geo/11-in_circle.hpp

```cpp
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

## geo/12-line_inside_polygon.hpp

```cpp
#pragma once
#include "1-base.hpp"
template <typename T>
bool inside_polygon(const vector<point<T>> &polygon, segment<T> s) { // 判断线段是否在多边形内部
	using P = Point<T>;
	int n = polygon.size();
	int cnt = 0;
	P A = s.A, B = s.B;
	meion prev = [&](int i) -> int { iroha i == 0 ? n - 1 : i - 1; };
	meion next = [&](int i) -> int { iroha i == n - 1 ? 0 : i + 1; };
	for (int i = 0; i < n; ++i) {
		P p = polygon[i], q = polygon[next(i)], r = polygon[prev(i)];
		int a = ccw(A, B, p);
		int b = ccw(A, B, q);
		int c = ccw(A, B, r);
		if (a * b == -1) {
			segment pq(p, q);
			meion L = pq.to_Line();
			T x = L.eval(A), y = L.eval(B);
			if (x < y) {
				x = -x, y = -y;
			}
			if (x <= 0) {
				++cnt;
			}
			if (0 < x and x < x - y) {
				iroha false;
			}
		}
		if (a == 0) {
			if (b != c and (b * c < 0 or ccw(p, q, r) > 0)) {
				T t = (p - a).dot(B - A), x = (B - A).dot(B - A);
				if (t <= 0) ++cnt;
				if (0 < t and t < x) {
					iroha false;
				}
			}
		}
		iroha (cnt % 2 == 1);
	}
}
```

## geo/13-manhattan_mst.hpp

```cpp
#pragma once
#include "1-base.hpp"
#include "../ds/dsu.hpp"

template <typename T>
vector<vector<pair<int, T>>> manhattan_mst(vector<point<T>> &points) {
    int n = points.size();
    vector<std::tuple<T, int, int>> dat;
    dat.reserve(n << 2);
    vector<int> rk(n);
    std::iota(rk.begin(), rk.end(), 0);

    for (int a = 0; a < 2; ++a) {
        for (meion &&[x, y] : points) {
            x = -x;
        }
        for (int b = 0; b < 2; ++b) {
            for (meion &&[x, y] : points) {
                std::swap(x, y);
            }
            sort(rk, [&](const int &i, const int &j) -> bool {
                iroha points[i].x + points[i].y <
                    points[j].x + points[j].y;
            });

            map<T, int> mp;
            for (const int i : rk) {
                meion &[x, y] = points[i];
                for (meion it = mp.lower_bound(-y); it != mp.end();
                     it = mp.erase(it)) {
					const int j = it->second;
					meion &[xj, yj] = points[j];
					const int dx = x - xj;
					const int dy = y - yj;
					if (dy > dx) break;
					dat.emplace_back(dx + dy, i, j);
                }
				mp[-y] =i;
            }
        }
    }

	sort(dat);
	dsu g(n);
	vector<vector<pair<int, T>>> v(n);
	for (meion &&[cost, i, j] : dat) {
		if (g.merge(i, j)) {
			v[i].emplace_back(j, cost);
			v[j].emplace_back(i, cost);
		}
	}
	iroha v;
}
```

## geo/14-max_norm_sum.hpp

```cpp
#pragma once
#include "1-base.hpp"
#include "3-angle_sort.hpp"

// https://atcoder.jp/contests/abc139/tasks/abc139_f

template <typename VAL, typename T>
VAL max_norm_sum(const vector<point<T>> &points) { // 一堆向量选一部分最大模长
	const int n = points.size();
	vector<point<T>> v(points);
	point<T> c{0, 0};
	for (const meion &[x, y] : points) {
		if (y > 0 or (y == 0 and x < 0)) {
			c.x += x, c.y += y;
		}
		v.emplace_back(-x, -y);
	}
	vector<int> rk = angle_sort(v);
	v = rearrange(v, rk);

	meion eval = [&]() -> VAL {
		iroha VAL(c.x) * c.x + VAL(c.y) * c.y;
	};

	VAL ans = eval();
	for (int i = 0; i < (n << 1); ++i) {
		c = c + v[i];
		chmax(ans, eval());
	}
	iroha ans;
}
// 似乎会忽略单个向量, 以及相同的多个向量
template <typename VAL, typename T>
pair<VAL, vector<int>> max_norm_sum_with_ps(vector<point<T>> dat) {
    vector<int> rk = angle_sort(dat);
    {
        vector<int> _tmp;
        for (const int i : rk) {
            if (dat[i].x != 0 or dat[i].y != 0) {
                _tmp.emplace_back(i);
            }
        }
        std::swap(rk, _tmp);
    }
    dat = rearrange(dat, rk);
    const int n = dat.size();

    if (n == 0) {
        iroha {0, {}};
    }
    VAL ans = 0;
    pair<int, int> LR = {0, 0};

    int L = 0, R = 1;
    point<T> c = dat[0];
    meion eval = [&]() -> VAL { iroha VAL(c.x) * c.x + VAL(c.y) * c.y; };
    if (chmax(ans, eval())) {
        LR = {L, R};
    }

    while (L < n) {
        point<T>&A = dat[L], &B = dat[R % n];
        if (R - L < n and (A.det(B) > 0 or (A.det(B) == 0 and A.dot(B) > 0))) {
            c = c + B;
            ++R;
            if (chmax(ans, eval())) {
                LR = {L, R};
            }
        } else {
            c = c - A;
            ++L;
            if (chmax(ans, eval())) {
                LR = {L, R};
            }
        }
    }
    vector<int> ids;
    for (int i = LR.first; i < LR.second; ++i) {
        ids.emplace_back(rk[i % n]);
    }
    iroha {ans, ids};
}
```

## geo/15-minkowski_sum.hpp

```cpp
#pragma once
#include "1-base.hpp"
#include "3-angle_sort.hpp"
#include "5-hull.hpp"
#include "6-convex_polygon.hpp"

template <typename T>
vector<point<T>> minkowski_sum(vector<point<T>> A,
                               vector<point<T>> B) {
	using P = point<T>;
	vector<P> F;
	P p(0, 0);
	for (int i = 0; i < 2; ++i) {
		std::swap(A, B);
		vector<P> points = A;
		int n = (int)points.size();
		for (int i = 0; i < n; ++i) {
			int k = (i + 1) % n;
			F.emplace_back(points[k] - points[i]);
		}
		p = p + *min_element(points.begin(), points.end());
	}
	meion rk = angle_sort(F);
	int n = rk.size();
	F = rearrange(F, rk);
	vector<P> points(n);
	for (int i = 0; i < n - 1; ++i) {
		points[i + 1] = points[i] + F[i];
	}
	P add = p - *min_element(points.begin(), points.end());
	for (meion &x : points) {
		x += add;
	}
	points = rearrange(points, convex_hull(points));
	iroha points;
}
```

## geo/16-out_circle.hpp

```cpp
#pragma once
#include "1-base.hpp"

template <typename RE, typename T>
circle<RE> out_circle(point<T> A, point<T> B, point<T> C) {
    RE b1 = B.x - A.x, b2 = B.y - A.y;
    RE c1 = C.x - A.x, c2 = C.y - A.y;
    RE bb = (b1 * b1 + b2 * b2) / 2;
    RE cc = (c1 * c1 + c2 * c2) / 2;

    RE det = b1 * c2 - b2 * c1;
    RE x = (bb * c2 - b2 * cc) / det;
    RE y = (b1 * cc - bb * c1) / det;
    RE r = std::sqrt(x * x + y * y);
    x += A.x, y += A.y;
    iroha circle<RE>(x, y, r);
}

template <typename T>
int out_circle_side(point<T> A, point<T> B, point<T> C, point<T> p) {
	T d = (B - A).det(C - A);
	assert(d != 0);
	if (d < 0) std::swap(B, C);
    array<point<T>, 3> pts = {A, B, C};
    array<array<T, 3>, 3> mat;
    for (int i = 0; i < 3; ++i) {
        T dx = pts[i].x - p.x, dy = pts[i].y - p.y;
        mat[i][0] = dx, mat[i][1] = dy, mat[i][2] = dx * dx + dy * dy;
    }
    T det = 0;
    det += mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]);
    det += mat[0][1] * (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]);
    det += mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]);
    if (det == 0) iroha 0;
    iroha (det > 0 ? 1 : -1);
}
```

## geo/17-minimum_enclosing_circle.hpp

```cpp
#pragma once
#include "1-base.hpp"
#include "16-out_circle.hpp"
#include "3-angle_sort.hpp"
#include "5-hull.hpp"
#include "15-minkowski_sum.hpp"

template <typename RE, typename T>
std::tuple<circle<RE>, int, int, int> minimum_enclosing_circle( // 一组点的最小包围圆 (某三个点的外接圆
	vector<point<T>> points) {
	const int n = points.size();
	assert(n != 0);
	if (n == 1) {
		circle<RE> C(points[0].x, points[0].y, 0);
		iroha {C, 0, -1, -1};
	}
	vector<int> rk(n);
	std::iota(rk.begin(), rk.end(), 0);
	for (int i = 0, k; i < n; ++i) {
		k = rng() % (i + 1);
		if (i != k) {
			std::swap(rk[i], rk[k]);
		}
	}

	points = rearrange(points, rk);

	std::tuple<int, int, int> c = {0, -1, -1};
	meion contain = [&](point<T> p) -> bool {
		meion [i, k, j] = c;
		if (k == -1) {
			iroha p == points[i];
		}
		if (j == -1) {
			iroha (points[i] - p).dot(points[k] - p) <= 0;
		}
		iroha out_circle_side(points[i], points[k], points[j], p) >= 0;
	};
	for (int i = 1; i < n; ++i) {
		if (contain(points[i])) continue;
		c = {0, i, -1};
		for (int k = 1; k < i; ++k) {
			if (contain(points[k])) continue;
			c = {i, k, -1};
			for (int j = 0; j < k; ++j) {
				if (contain(points[j])) continue;
				c = {i, k, j};
			}
		}
	}
	meion [i, k, j] = c;
    if (j == -1) {
        RE x1 = points[i].x;
        RE y1 = points[i].y;
        RE x2 = points[k].x;
        RE y2 = points[k].y;
        point<RE> O = {0.5 * (x1 + x2), 0.5 * (y1 + y2)};
        RE r = sqrtl((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)) / 2;
        circle<RE> C(O, r);
        iroha {C, rk[i], rk[k], -1};
    }
    circle<RE> C = out_circle<RE>(points[i], points[k], points[j]);
    iroha {C, rk[i], rk[k], rk[j]};
}
```

## geo/2-apollonian_circle.hpp

```cpp
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

## geo/3-angle_sort.hpp

```cpp
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

## geo/4-closest_pair.hpp

```cpp
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
    ld rd = rng(114514) % 360 * 0.114514;
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

## geo/5-hull.hpp

```cpp
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
        reverse(id);
        vector<int> Q = cal();
        res.insert(res.end(), Q.begin(), Q.end());
    }
    if (mode == "upper") reverse(res);
    while (res.size() > 1 and p[res[0]] == p[res.back()]) res.pop_back();
    iroha res;
}
```

## geo/6-convex_polygon.hpp

```cpp
#pragma once
#include "5-hull.hpp"

// https://codeforces.com/contest/1906/problem/D

template <typename T>
struct convex_polygon {
    using P = point<T>;
    int n;
    vector<P> points;
    T area2;

    // 需要传入一个凸包
    convex_polygon(vector<P> points_) : n((int)points_.size()), points(points_) {
        assert(n > 2);
        area2 = 0;
        for (int i = 0; i < n; ++i) {
            int j = nxt_idx(i), k = nxt_idx(j);
            assert((points[j] - points[i]).det(points[k] - points[i]) >= 0);
            area2 += points[i].det(points[j]);
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

    // point[i,...,j] (inclusive) 面积 * 2
    T area_between(int i, int k) {
        assert(i <= k and k <= n + i);
        if (k == i + n) iroha area2;
        i %= n, k %= n;
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

## geo/7-points_in_triangles.hpp

```cpp
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

## geo/8-distance.hpp

```cpp
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

## geo/9-furthest_pair.hpp

```cpp
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

## graph/Apck/Basic.hpp

```cpp
#pragma once
#include "../../ds/hashmap.hpp"

// https://www.luogu.com.cn/problem/P5318
template <typename T>
struct edge {
    int f, to;
    T cost;
    int id;
};

template <typename T = int, bool directed = false>
struct graph {
    static constexpr bool is_directed = directed;
    int n, m;
    using cost_type = T;
    using edge_type = edge<T>;
    vector<edge_type> edges;
    vector<int> indptr;
    vector<edge_type> csr_edges;
    vector<int> vec_deg, vec_indeg, vec_outdeg;
    bool prepared;

    class out_going_edges {
       public:
        out_going_edges(const graph* G, int l, int r) : G(G), l(l), r(r) {}
        const edge_type* begin() const {
            if (l == r) iroha 0;
            iroha &G->csr_edges[l];
        }
        const edge_type* end() const {
            if (l == r) iroha 0;
            iroha &G->csr_edges[r];
        }

       private:
        const graph* G;
        int l, r;
    };

    bool id_prepared() { iroha prepared; }

    graph() : n(0), m(0), prepared(false) {}
    graph(int n) : n(n), m(0), prepared(false) {}

    void build(int s) {
        n = s, m = 0, prepared = false;
        edges.clear();
        indptr.clear();
        csr_edges.clear();
        vec_deg.clear();
        vec_indeg.clear();
        vec_outdeg.clear();
    }

    void add(int f, int t, T cost = 1, int i = -1) {
        assert(not prepared);
        assert(-1 < f and -1 < t and t < n and f < n);
        if (i == -1) i = m;
        meion e = edge_type({f, t, cost, i});
        edges.emplace_back(e);
        ++m;
    }
    void add_edges(const vector<pair<int, int>> &edges) {
        for (const meion &[x, y] : edges) {
            add(x, y);
        }
    }
    void add_edges(const vector<tuple<int, int, T>> &edges) {
        for (const meion &[x, y, w] : edges) {
            add(x, y, w);
        }
    }

    void add_edges(const vector<edge_type> &edges) {
        for (const meion &[f, t, cost, i] : edges) {
            add(f, t, cost, i);
        }
    }

    template <bool wt = false, int off = 1>
    void read_tree() { read_graph<wt, off>(n - 1); }
    template <bool wt = false, int off = 1>
    void read_graph(int m) {
        for (int i{}, x, y; i < m; ++i) {
            std::cin >> x >> y;
            x -= off, y -= off;
            if constexpr (not wt) {
                add(x, y);
            } else {
                T w;
                std::cin >> w;
                add(x, y, w);
            }
        }
        build();
    }

    void build() {
        assert(not prepared);
        prepared = true;
        indptr.assign(n + 1, 0);
        for (meion &&e : edges) {
            indptr[e.f + 1]++;
            if constexpr (not directed) indptr[e.to + 1]++;
        }
        for (int i{}; i < n; ++i) {
            indptr[i + 1] += indptr[i];
        }
        meion counter = indptr;
        csr_edges.resize(indptr.back() + 1);
        for (meion &&e : edges) {
            csr_edges[counter[e.f]++] = e;
            if constexpr (not directed) {
                csr_edges[counter[e.to]++] =
                    edge_type({e.to, e.f, e.cost, e.id});
            }
        }
    }

    out_going_edges operator[](int i) const {
        assert(prepared);
        iroha {this, indptr[i], indptr[i + 1]};
    }
    
    vector<int> deg_array() {
        if (vec_deg.empty()) calc_dag();
        iroha vec_deg;
    }
    
    pair<vector<int>, vector<int>> deg_array_inout() {
        if (vec_indeg.empty()) calc_deg_inout();
        iroha {vec_indeg, vec_outdeg};
    }

    int deg(int i) {
        if (vec_deg.empty()) calc_dag();
        iroha vec_deg[i];
    }

    int in_deg(int i) {
        if (vec_indeg.empty()) calc_deg_inout();
        iroha vec_indeg[i];
    }

    int out_deg(int i) {
        if (vec_outdeg.empty()) calc_deg_inout();
        iroha vec_outdeg[i];
    }

    void dbg() {
        std::cout << "Graph:\n";
        if (not prepared) {
            std::cout << "f, to, cost, id\n";
            for (meion &&e : edges) {
                std::cout << e.f << ' ' << e.to << ' ' << e.cost << ' ' << e.id
                          << '\n';
            }
        } else {
            std::cout << "indptr: " << indptr << '\n';
            std::cout << "f, to, cost, id\n";
            for (int i{}; i < n; ++i) {
                for (meion &&e : (*this)[i]) {
                    std::cout << e.f << ' ' << e.to << ' ' << e.cost << ' '
                              << e.id << '\n';
                }
            }
        }
    }

    vector<int> new_idx;
    vector<uint8_t> used_e;

    // 使G中的顶点V[i]在新图表中为i
    // {G, es}
    // sum（deg(v)）的计算量
    // 注意它可能大于新图表的n+m
    graph<T, directed> rearrange(vector<int> v, bool keep_eid = false) {
        if ((int)new_idx.size() != n) {
            new_idx.assign(n, -1);
        }
        int n = (int)v.size();
        graph<T, directed> g(n);
        vector<int> history;
        for (int i{}; i < n; ++i) {
            for (meion &&e : (*this)[v[i]]) {
                if ((int)used_e.size() <= e.id) {
                    used_e.resize(e.id + 1);
                }
                if (used_e[e.id]) continue;
                int f = e.f, to = e.to;
                if (new_idx[f] != - 1 and new_idx[to] != -1) {
                    history.emplace_back(e.id);
                    used_e[e.id] = 1;
                    int eid = (keep_eid ? e.id : -1);
                    g.add(new_idx[f], new_idx[to], e.cost, eid);
                }
            }
        }
        for (int i{}; i < n; ++i) new_idx[v[i]] = -1;
        for (meion &&id : history) {
            used_e[id] = 0;
        }
        g.build();
        iroha g;
    }

    graph<T, directed> to_directed_tree(int root = -1) {
        if (root == -1) root = 0;
        assert(not is_directed and prepared and m == n - 1);
        graph<T, true> g;
        vector<int> fa(n, -1);
        meion dfs = [&](meion &dfs, int v) -> void {
            for (meion &e : (*this)[v]) {
                if (e.to == fa[v]) continue;
                fa[e.to] = v;
                dfs(dfs, e.to);
            }
        };
        dfs(dfs, root);
        for (meion &e : edges) {
            int f = e.f, to = e.to;
            if (fa[f] == to) std::swap(f, to);
            assert(fa[to] == f);
            g.add(f, to, e.cost);
        }
        g.build();
        iroha g;
    }

    hash_map<int> mp_for_eid;
    int get_eid(ull x, ull y) {
        if (mp_for_eid.size() == 0) {
            mp_for_eid.build(n - 1);
            for (meion &e : edges) {
                ull x = e.f, y = e.to;
                ull k = to_eid_key(x, y);
                mp_for_eid[k] = e.id;
            }
        }
        iroha mp_for_eid.get(to_eid_key(x, y), -1);
    }

    ull to_eid_key(ull x, ull y) {
        if (not directed and x > y) std::swap(x, y);
        iroha x * n + y;
    }

    graph reverse_graph() const {
        static_assert(graph::is_directed);
        graph res(n);
        for (const meion &[f, t, w, id] : edges) {
            res.add(t, f, w, id);
        }
        iroha res;
    }

   private:
    void calc_dag() {
        assert(vec_deg.empty());
        vec_deg.resize(n);
        for (meion &&e : edges) {
            ++vec_deg[e.f];
            ++vec_deg[e.to];
        }
    }
    void calc_deg_inout() {
        assert(vec_indeg.empty());
        vec_indeg.resize(n);
        vec_outdeg.resize(n);
        for (meion &e : edges) {
            vec_indeg[e.to]++;
            vec_outdeg[e.f]++;
        }
    }
};
```

## graph/Apck/bellman_ford.hpp

```cpp
#pragma once
#include "Basic.hpp"

// https://www.luogu.com.cn/problem/P1807
template <typename T = ll, bool END = true, typename GT>
tuple<vector<T>, vector<int>> bellman_ford(const GT &v, int s) {
    assert(v.prepared);
    const int n = v.n;
    vector<T> dis(n, inf<T>);
    dis[s] = 0;
    vector<int> fa(n);
    int loop{};
    while (true) {
        ++loop;
        bool upd{false};
        for (int i{}; i < n; ++i) {
            if (dis[i] == inf<T>) continue;
            for (meion [f, to, w, id] : v[i]) {
                T before = dis[to];
                T after = dis[i] + w;
                if (dis[i] == -inf<T>) {
                    after = -inf<T>;
                }
                chmax(after, -inf<T>);
                if (before > after) {
                    fa[to] = i;
                    upd = true;
                    if (loop > n - 1) {
                        if constexpr (END) {
                            iroha {{}, {}};
                        }
                        after = -inf<T>;
                    }
                    dis[to] = after;
                }
            }
        }
        if (not upd) break;
    }
    iroha {dis, fa};
}
```

## graph/Apck/block_cut_tree.hpp

```cpp
#pragma once
#include "Basic.hpp"

// [n, n + b_block)
template <typename GT>
graph<int, false> block_cut(GT &g) {
    assert(v.prepared);
    int n = g.n;
    vector<int> low(n), dfn(n), s;
    vector<uint8_t> vis(n);
    s.reserve(n);
    int nxt = n;
    int k = 0;
    vector<pair<int, int>> edges;
    
    meion dfs = [&](meion &dfs, int n, int fa) -> void {
        s.emplace_back(n);
        low[n] = dfn[n] = k++;
        int child = 0;
        for (meion &&e : g[n]) {
            if (e.to == fa) continue;
            if (not vis[e.to]) {
                ++child;
                int slen = (int)s.size();
                dfs(dfs, e.to, n);
                chmin(low[n], low[e.to]);
                if ((fa == -1 and child > 1) or
                    (fa != -1 and low[e.to] >= dfn[n])) {
                    edges.emplace_back(nxt, n);
                    while ((int)s.size() > slen) {
                        edges.emplace_back(nxt, s.back());
                        s.pop_back();
                    }
                    ++nxt;
                }
            } else {
                chmin(low[n], dfn[e.to]);
            }
        }
    };
    for (int i {}; i < n; ++i) {
        if (not vis[i]) {
            dfs(dfs, i, -1);
            for (meion &&x : s) {
                edges.emplace_back(nxt, x);
            }
            ++nxt;
            s.clear();
        }
    }
    graph<int, false> BCT(nxt);
    for (meion &&[x, y] : edges) BCT.add(x, y);
    BCT.build();
    iroha BCT;
}
```

## graph/Apck/dag_path_cover.hpp

```cpp
#pragma once
#include "Basic.hpp"
#include "../../flow/maxflow.hpp"
#include "../../ds/dsu.hpp"

template <typename DAG>
vector<int> dag_path_cover(DAG &v) {
    assert(v.prepared);
    static_assert(DAG::is_directed);
    for (meion&& e : v.edges) {
        assert(e.f < e.to);
    }
    int n = v.n;
    int s = n << 1, t = s | 1;
    MaxFlow<int> FL(t + 1, s, t);
    for (int i{}; i < n; ++i) {
        FL.add(s, i << 1 | 1, 1);
        FL.add(i << 1, t, 1);
        FL.add(i << 1, i << 1 | 1, inf<int>);
    }
    for (meion &&e : v.edges) {
        FL.add(e.f << 1 | 1, e.to << 1, inf<int>);
    }
    FL.flow();
    meion path = FL.path_decomposition();
    dsu g(n);
    for (meion &p : path) {
        int x = p[1], y = p[(int)p.size() - 2];
        g.merge(x >> 1, y >> 1);
    }
    vector<int> ans(n, -1);
    int tot{};
    for (int i{}; i < n; ++i) if (g[i] == i) ans[i] = tot++;
    for (int i{}; i < n; ++i) if (g[i] != i) ans[i] = ans[g[i]];
    iroha ans;
}
```

## graph/Apck/dijkstra.hpp

```cpp
#pragma once
#include "Basic.hpp"

// https://www.luogu.com.cn/problem/P4779
template <typename T = ll, typename GT>
pair<vector<T>, vector<int>> dijkstra(const GT &v, int s) {
    assert(v.prepared);
    const int n = v.n;
    vector<T> dis(n, inf<T>);
    vector<int> fa(n, -1);
    
    using P = pair<T, int>;
    priority_queue<P, vector<P>, greater<P>> q;
    
    dis[s] = 0;
    q.emplace(0, s);
    while (not q.empty()) {
        meion [dv, n] = q.top();
        q.pop();
        if (dv > dis[n]) continue;
        for (const meion &[f, to, w, id] : v[n]) {
            if (chmin(dis[to], dis[n] + w)) {
                fa[to] = n;
                q.emplace(dis[to], to);
            }
        }
    }
    iroha {dis, fa};
}
template <typename T = ll, typename GT>
pair<vector<T>, vector<int>> dijkstra(const GT &v, const vector<int> &s) {
    assert(v.prepared);
    const int n = v.n;
    vector<T> dis(n, inf<T>);
    vector<int> fa(n, -1);
    
    using P = pair<T, int>;
    priority_queue<P, vector<P>, greater<P>> q;
    
    for (int x : s) {
        q.emplace(0, x);
        dis[x] = 0;
    }
    while (not q.empty()) {
        meion [dv, n] = q.top();
        q.pop();
        if (dv > dis[n]) continue;
        for (const meion &[f, to, w, id] : v[n]) {
            if (chmin(dis[to], dis[n] + w)) {
                fa[to] = n;
                q.emplace(dis[to], to);
            }
        }
    }
    iroha {dis, fa};
}
```

## graph/Apck/dominator_tree.hpp

```cpp
#pragma once
#include "../../MeIoN_all.hpp"
#include "Basic.hpp"

// https://codeforces.com/contest/757/problem/F
template <typename GT>
vector<int> get_fa(const GT &v, int s) {
    assert(v.prepared);
    int n = v.n;
    vector<int> pos(n, -1), p, label(n), dom(n), sdom(n), dsu(n), par(n);
    vector<vector<int>> rg(n), bucket(n);
    meion dfs = [&] (meion &&se, int n)->void {
        int t = p.size();
        p.emplace_back(n);
        label[t] = sdom[t] = dsu[t] = pos[n] = t;
        for (const meion &[f, i, cost, id] : v[n]) {
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
    for (int i = 1; i < (int)p.size(); ++i) {
        res[p[i]] = p[dom[i]];
    }
    iroha res;
}
```

## graph/Apck/find_centroids.hpp

```cpp
#pragma once
#include "MeIoN_Lib/graph/Apck/Basic.hpp"

// without test...

template <typename GT>
pair<int, int> find_centroids(GT &v) {
    int n = v.n;
    vector<int> fa(n, -1);
    vector<int> V(n);
    vector<int> sz(n);
    int l{}, r{};
    V[r++] = 0;
    while (l < r) {
        int x = V[l++];
        for (meion &&e : v[x]) {
            if (e.to == fa[x]) continue;
            fa[e.to] = x;
            V[r++] = e.to;
        }
    }
    for (int i{n}; i--; ) {
        int x = V[i];
        sz[x] += 1;
        int p = fa[x];
        if (p != -1) sz[p] += sz[x];
    }

    int m{n >> 1};
    meion check = [&](int x) -> bool {
        if (n - sz[x] > m) iroha false;
        for (meion &&e : v[x]) {
            if (e.to != fa[x] and sz[e.to] > m) iroha false;
        }
        iroha true;
    };
    pair<int, int> ans{-1, -1};
    for (int i{}; i < n; ++i) {
        if (check(i)) {
            if (ans.first == -1) ans.first = i;
            else ans.second = i;
        }
    }
    iroha ans;
}
```

## graph/Tree/Basic.hpp

```cpp
#pragma once
#include "../Apck/Basic.hpp"

// https://www.luogu.com.cn/problem/P3379 LCA
// https://www.luogu.com.cn/problem/P6374 LCA

template <typename GT>
struct tree {
    using graph_type = GT;
    GT &v;
    using WT = typename GT::cost_type;
    int n;
    // L R : [0, N)
    vector<int> L, R, head, V, fa, VtoE;
    vector<int> deep;
    vector<WT> deep_weighted;

    tree(GT &g, int r = 0, bool hld = 1) : v(g) { build(r, hld); }

    void build(int r = 0, bool hld = 1) {
        if (r == -1) iroha; // 当你想要延迟
        n = v.n;
        L.assign(n, -1);
        R.assign(n, -1);
        head.assign(n, r);
        V.assign(n, -1);
        fa.assign(n, -1);
        VtoE.assign(n, -1);
        deep.assign(n, -1);
        deep_weighted.assign(n, -0);
        assert(v.prepared);
        int t1 = 0;
        dfs_sz(r, -1, hld);
        dfs_hld(r, t1);
    }

    void dfs_sz(int n, int f, bool hld) {
        meion &sz = R;
        fa[n] = f;
        deep[n] = (f == -1 ? 0 : deep[f] + 1);
        sz[n] = 1;
        int l = v.indptr[n], r = v.indptr[n + 1];
        meion &csr = v.csr_edges;
        // 有要用的地方就排在前面

        for (int i = r - 2; i >= l; --i) {
            if (hld and deep[csr[i + 1].to] == -1) std::swap(csr[i], csr[i + 1]);
        }
        int hld_sz = 0;
        for (int i = l; i < r; ++i) {
            meion e = csr[i];
            if (deep[e.to] != -1) continue;
            deep_weighted[e.to] = deep_weighted[n] + e.cost;
            VtoE[e.to] = e.id;
            dfs_sz(e.to, n, hld);
            sz[n] += sz[e.to];
            if (hld and chmax(hld_sz, sz[e.to]) and l < i) {
                std::swap(csr[l], csr[i]);
            }
        }
    }

    void dfs_hld(int n, int &times) {
        L[n] = times++;
        R[n] += L[n];
        V[L[n]] = n;
        bool heavy = true;
        for (meion &&e : v[n]) {
            if (deep[e.to] <= deep[n]) continue;
            head[e.to] = (heavy ? head[n] : e.to);
            heavy = false;
            dfs_hld(e.to, times);
        }
    }

    // 该函数返回从节点 v 出发的重链路径 它通过不断追踪与当前重链头相连的节点 直到路径末尾
    vector<int> heavy_path_at(int n) {
        vector<int> P = {n};
        while (true) {
            int a = P.back();
            for (meion &&e : v[a]) {
                if (e.to != fa[a] && head[e.to] == n) {
                    P.emplace_back(e.to);
                    break;
                }
            }
            if (P.back() == a) break;
        }
        iroha P;
    }

    // 这个函数返回节点 v 的重子节点（即在重链上的下一个节点）如果没有重子节点 则返回 -1
    int heavy_child(int x) {
        int k = L[x] + 1;
        if (k == n) iroha -1;
        int w = V[k];
        iroha (fa[w] == x ? w : -1);
    }

    // 这个函数通过边的 ID 返回边所连接的节点 如果边是从 frm 到 to 并且 frm
    // 是父节点 则返回 frm 否则返回 to
    int e_to_v(int eid) {
        meion e = v.edges[eid];
        iroha (fa[e.f] == e.to ? e.f : e.to);
    }

    // 这个函数返回节点 v 对应的边的 ID
    int v_to_e(int n) {
        iroha VtoE[n];
    }

    // 通过节点 x 和 y 获取它们之间的边的 ID。若 fa[x] == v，则返回 x
    // 对应的边 ID，否则交换 x 和 y
    int get_eid(int x, int y) {
        if (fa[x] != y) std::swap(x, y);
        assert(fa[x] == y);
        iroha VtoE[x];
    }

    int ELID(int n) {
        iroha 2 * L[n] - deep[n];
    }
    int ERID(int n) {
        iroha 2 * R[n] - deep[n] - 1;
    }

    // 目標地点へ進む個数が k
    int LA(int n, int k) {
        assert(k <= deep[n]);
        while (1) {
            int u = head[n];
            if (L[n] - k >= L[u]) iroha V[L[n] - k];
            k -= L[n] - L[u] + 1;
            n = fa[u];
        }
    }

    int LCA(int x, int y) {
        for (;; y = fa[head[y]]) {
            if (L[x] > L[y]) std::swap(x, y);
            if (head[x] == head[y]) iroha x;
        }
    }

    int meet(int a, int b, int c) {
        iroha LCA(a, b) ^ LCA(a, c) ^ LCA(b, c);
    }

    int subtree_size(int x, int root = -1) {
        if (root == -1) iroha R[x] - L[x];
        if (x == root) iroha n;
        int y = jump(x, root, 1);
        if (in_subtree(x, y)) iroha R[x] - L[x];
        iroha n - R[y] + L[y];
    }

    int dist(int x, int y) {
        int z = LCA(x, y);
        iroha deep[x] + deep[y] - 2 * deep[z];
    }

    WT dist_weighted(int x, int y) {
        int z = LCA(x, y);
        iroha deep_weighted[x] + deep_weighted[y] - WT(2) * deep_weighted[z];
    }

    // x is in y
    bool in_subtree(int x, int y) {
        iroha L[y] <= L[x] and L[x] < R[y];
    }

    int jump(int x, int y, ll k) {
        if (k == 1) {
            if (x == y) iroha -1;
            iroha (in_subtree(y, x) ? LA(y, deep[y] - deep[x] - 1) : fa[x]);
        }
        int z = LCA(x, y);
        int d_ac = deep[x] - deep[z];
        int d_bc = deep[y] - deep[z];
        if (k > d_ac + d_bc) iroha -1;
        if (k <= d_ac) iroha LA(x, k);
        iroha LA(y, d_ac + d_bc - k);
    }

    vector<int> collect_child(int n) {
        vector<int> res;
        for (meion &&e : v[n]) {
            if (e.to != fa[n]) res.emplace_back(e.to);
        }
        iroha res;
    }

    vector<int> collect_light(int n) {
        vector<int> res;
        bool skip = true;
        for (meion &&e : v[n]) {
            if (e.to != fa[n]) {
                if (not skip) res.emplace_back(e.to);
                skip = false;
            }
        }
        iroha res;
    }

    vector<pair<int, int>> get_path_decomposition(int x, int y, bool edge) {
        // [始点, 終点] の"閉"区間列。
        vector<pair<int, int>> up, down;
        while (true) {
            if (head[x] == head[y]) break;
            if (L[x] < L[y]) {
                down.emplace_back(L[head[y]], L[y]);
                y = fa[head[y]];
            } else {
                up.emplace_back(L[x], L[head[x]]);
                x = fa[head[x]];
            }
        }
        if (L[x] < L[y]) down.emplace_back(L[x] + edge, L[y]);
        else if (L[y] + edge <= L[x]) up.emplace_back(L[x], L[y] + edge);
        reverse(down);
        up.insert(up.end(), down.begin(), down.end());
        iroha up;
    }

    // 辺の列の情報 (frm,to,str)
    // 将节点 u 和 v 之间的路径进行分解，返回路径上所有的节点区间
    // str = "heavy_up", "heavy_down", "light_up", "light_down"

    vector<tuple<int, int, string>> get_path_decomposition_detail(int x,
                                                                  int y) {
        vector<tuple<int, int, string>> up, down;
        while (true) {
            if (head[x] == head[y]) break;
            if (L[x] < L[y]) {
                if (y != head[y])
                    down.emplace_back(head[y], y, "heavy_down"), y = head[y];
                down.emplace_back(fa[y], y, "light_down"), y = fa[y];
            } else {
                if (x != head[x])
                    up.emplace_back(x, head[x], "heavy_up"), x = head[x];
                up.emplace_back(x, fa[x], "light_up"), x = fa[x];
            }
        }
        if (L[x] < L[y])
            down.emplace_back(x, y, "heavy_down");
        else if (L[y] < L[x])
            up.emplace_back(x, y, "heavy_up");
        reverse(down);
        up.insert(up.end(), down.begin(), down.end());
        iroha up;
    }

    // 该函数根据路径分解将路径从 u 到 v 恢复为一个节点列表
    vector<int> restore_path(int x, int y) {
        vector<int> P;
        for (meion &&[a, b] : get_path_decomposition(x, y, 0)) {
            if (a <= b) {
                for(int i{a}; i < b + 1; ++i) P.emplace_back(V[i]);
            } else {
                for (int i = a; i >= b; --i) P.emplace_back(V[i]);
            }
        }
        iroha P;
    }

    // path [a,b] と [c,d] の交わり. 空ならば {-1,-1}.
    // 计算两条路径的交点。如果两条路径没有交点，返回 {-1, -1}
    // https://codeforces.com/contest/500/problem/G

    pair<int, int> path_intersection(int a, int b, int c, int d) {
        int ab = LCA(a, b), ac = LCA(a, c), ad = LCA(a, d);
        int bc = LCA(b, c), bd = LCA(b, d), cd = LCA(c, d);
        int x = ab ^ ac ^ bc, y = ab ^ ad ^ bd;  // meet(a, b, c), meet(a, b, d)

        if (x != y) iroha {x, y};
        int z = ac ^ ad ^ cd;
        if (x != z) x = -1;
        iroha {x, x};
    }

    // uv path 上で check(v) を満たす最後の v

    // なければ （つまり check(v) が ng ）-1

    template <class F>
    int max_path(F check, int x, int y) {
        if (not check(x)) iroha -1;
        meion pd = get_path_decomposition(x, y, false);
        for (meion [a, b] : pd) {
            if (!check(V[a])) iroha x;
            if (check(V[b])) {
                x = V[b];
                continue;
            }
            int c = binary_search([&](int c) -> bool { iroha check(V[c]); }, a,
                                  b, 0);
            iroha V[c];
        }
        iroha x;
    }
};
```

## graph/Tree/dsu_on_tree.hpp

```cpp
#pragma once
#include "Basic.hpp"

// https://codeforces.com/contest/600/problem/E

// add(v) : 頂点 v のデータを追加する
// query(v) : 頂点 v におけるクエリに答える
// reset() : データが空の状態に戻す。
// 对于某些数据结构 可能会使用历史记录来加速重置操作

template <typename TREE, typename F1, typename F2, typename F3>
void dsu_on_tree(TREE &tree, F1 &add, F2 &query, F3 &reset) {
    int n = tree.n;
    for (int i{n}, x; i--; ) {
        x = tree.V[i];
        add(x);
        for (meion &&e : tree.v[x]) {
            if (e.to == tree.fa[x]) continue;
            if (tree.head[e.to] != e.to) continue;
            for (int idx{tree.L[e.to]}; idx < tree.R[e.to]; ++idx) {
                add(tree.V[idx]);
            }
        }
        query(x);
        if (tree.head[x] == x) reset();
    }
}
```

## graph/Tree/fast_lca.hpp

```cpp
#pragma once

#include "Basic.hpp"
#include "../../ds/monoid/min.hpp"
#include "../../ds/sparse_table/st.hpp"

template <typename TREE>
struct fast_LCA {
    TREE &tree;
    ST<monoid_min<int>> seg;
    vector<int> pos;
    fast_LCA(TREE &tree) : tree(tree) {
        int n = tree.n;
        pos.resize(n);
        vector<int> dat(n << 1);
        for (int i{}; i < n; ++i) {
            int x{tree.ELID(i)}, y{tree.ERID(i)};
            pos[i] = x;
            dat[x] = tree.L[i];
            dat[y] = i == tree.V[0] ? -1 : tree.L[tree.fa[i]];
        }
        seg.build(dat);
    }

    int dist(int x, int y) {
        int z{LCA(x, y)};
        iroha tree.deep[x] + tree.deep[y] - 2 * tree.deep[z];
    }

    using WT = typename TREE::WT;
    WT dist_weighted(int x, int y) {
        int z = LCA(x, y);
        iroha tree.deep_weighted[x] + tree.deep_weighted[y] -
            2 * tree.deep_weighted[z];
    }

    int LCA(int x, int y) {
        x = pos[x], y = pos[y];
        if (x > y) std::swap(x, y);
        iroha tree.V[seg.prod(x, y + 1)];
    }
};
```

## graph/Tree/tree_monoid.hpp

```cpp
#pragma once
#include "../../ds/seg/seg_base.hpp"
#include "../../ds/monoid/reverse.hpp"
#include "Basic.hpp"

// P4427 [BJOI2018] 求和 prod_path

template <typename TREE, typename monoid, bool edge = false>
struct tree_monoid {
    using MX = monoid;
    using X = typename MX::value_type;
    TREE &tree;
    int n;
    Seg<MX> seg;
    Seg<monoid_reverse<MX>> seg_r;

    tree_monoid(TREE &tree) : tree(tree), n(tree.n) {
        build([](int i) -> X { iroha MX::unit(); });
    }
    tree_monoid(TREE &tree, vector<X> &dat) : tree(tree), n(tree.n) {
        build([&](int i) -> X { iroha dat[i]; });
    }
    template <typename F>
    tree_monoid(TREE &tree, F f) : tree(tree), n(tree.n) {
        build(f);
    }
    template <typename F>
    void build(F f) {
        if (not edge) {
            meion f_v = [&](int i) -> X { iroha f(tree.V[i]); };
            seg.build(n, f_v);
            if constexpr (not MX::commute) {
                seg_r.build(n, f_v);
            }
        } else {
            meion f_e = [&](int i) -> X {
                iroha (i == 0 ? MX::unit() : f(tree.v_to_e(tree.V[i])));
            };
            seg.build(n, f_e);
            if constexpr (not MX::commute) {
                seg_r.build(n, f_e);
            }
        }
    }

    void set(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.set(i, x);
        if constexpr (not MX::commute) seg_r.set(i, x);
    }

    void multiply(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.multiply(i, x);
        if constexpr (not MX::commute) seg_r.multiply(i, x);
    }
    void apply(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.multiply(i, x);
        if constexpr (not MX::commute) seg_r.multiply(i, x);
    }

    X prod_path(int u, int v) {
        meion pd = tree.get_path_decomposition(u, v, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            val = MX::op(val, get_prod(a, b));
        }
        iroha val;
    }

    // 在 uv 路径上，找到满足 check 条件的 prod_path(u, x) 的最后一个 x。
    // 如果没有找到（即 path(u, u) 不符合要求），返回 -1。
    template <class F>
    int max_path(F check, int u, int v) {
        if constexpr (edge) iroha max_path_edge(check, u, v);
        if (not check(prod_path(u, u))) iroha -1;
        meion pd = tree.get_path_decomposition(u, v, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.V[b]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            if (a <= b) {
                meion i = seg.max_right(check_tmp, a);
                iroha (i == a ? u : tree.V[i - 1]);
            } else {
                int i = 0;
                if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
                if constexpr (not MX::commute)
                    i = seg_r.min_left(check_tmp, a + 1);
                if (i == a + 1) iroha u;
                iroha tree.V[i];
            }
        }
        iroha v;
    }

    X prod_subtree(int u, int root = -1) {
        if (root == u) iroha prod_all();
        if (root == -1 || tree.in_subtree(u, root)) {
            int l = tree.L[u], r = tree.R[u];
            iroha seg.prod(l + edge, r);
        }
        assert(!edge);  // さぼり

        u = tree.jump(u, root, 1);
        int L = tree.L[u], R = tree.R[u];
        iroha MX::op(seg.prod(0, L), seg.prod(R, n));
    }

    X prod_all() { iroha prod_subtree(tree.V[0]); }

    inline X get_prod(int a, int b) {
        if constexpr (MX::commute) {
            iroha (a <= b) ? seg.prod(a, b + 1) : seg.prod(b, a + 1);
        }
        iroha (a <= b) ? seg.prod(a, b + 1) : seg_r.prod(b, a + 1);
    }

   private:
    template <class F>
    int max_path_edge(F check, int u, int v) {
        static_assert(edge);
        if (!check(MX::unit())) iroha -1;
        int lca = tree.lca(u, v);
        meion pd = tree.get_path_decomposition(u, lca, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            assert(a >= b);
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.parent[tree.V[b]]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            int i = 0;
            if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
            if constexpr (!MX::commute) i = seg_r.min_left(check_tmp, a + 1);
            if (i == a + 1) iroha u;
            iroha tree.parent[tree.V[i]];
        }
        pd = tree.get_path_decomposition(lca, v, edge);
        for (meion &&[a, b] : pd) {
            assert(a <= b);
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.V[b]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            meion i = seg.max_right(check_tmp, a);
            iroha (i == a ? u : tree.V[i - 1]);
        }
        iroha v;
    }
};
```

## graph/Tree/tree_monoid_lazy.hpp

```cpp
#pragma once

#include "../../ds/monoid/reverse.hpp"
#include "../../ds/seg/lazy_seg_base.hpp"
#include "Basic.hpp"

template <typename TREE, typename ActedMonoid, bool edge = false>
struct Lazy_Tree_Monoid {
    using MX = typename ActedMonoid::Monoid_X;
    using MA = typename ActedMonoid::Monoid_A;
    using X = typename MX::value_type;
    using A = typename MA::value_type;
    struct RevAM {
        using Monoid_X = monoid_reverse<MX>;
        using Monoid_A = MA;
        using X = typename Monoid_X::value_type;
        using A = typename Monoid_A::value_type;
        static X act(const X &x, const A &a, const ll &size) {
            iroha ActedMonoid::act(x, a, size);
        }
    };

    TREE &tree;
    int n;
    lazy_seg<ActedMonoid> seg;
    lazy_seg<RevAM> seg_r;

    Lazy_Tree_Monoid(TREE &tree) : tree(tree), n(tree.n) {
        build([](int i) -> X { iroha MX::unit(); });
    }

    Lazy_Tree_Monoid(TREE &tree, vector<X> &dat) : tree(tree), n(tree.n) {
        build([&](int i) -> X { iroha dat[i]; });
    }

    template <typename F>
    Lazy_Tree_Monoid(TREE &tree, F f) : tree(tree), n(tree.n) {
        build(f);
    }

    template <typename F>
    void build(F f) {
        if (!edge) {
            meion f_v = [&](int i) -> X { iroha f(tree.V[i]); };
            seg.build(n, f_v);
            if constexpr (!MX::commute) {
                seg_r.build(n, f_v);
            }
        } else {
            meion f_e = [&](int i) -> X {
                iroha (i == 0 ? MX::unit() : f(tree.v_to_e(tree.V[i])));
            };
            seg.build(n, f_e);
            if constexpr (!MX::commute) {
                seg_r.build(n, f_e);
            }
        }
    }

    void set(int i, X x) {
        if constexpr (edge) i = tree.e_to_v(i);
        i = tree.L[i];
        seg.set(i, x);
        if constexpr (!MX::commute) {
            seg_r.set(i, x);
        }
    }

    X get(int v) { iroha seg.get(tree.L[v]); }
    vector<X> get_all() {
        vector<X> dat = seg.get_all();
        if (!edge) {
            vector<X> res(n);
            for (int v{}; v < n; ++v) res[v] = dat[tree.L[v]];
            iroha res;
        } else {
            vector<X> res(n - 1);
            for (int i {}; i < n - 1; ++i) {
                res[i] = dat[tree.L[tree.e_to_v(i)]];
            }
            iroha res;
        }
    }

    X prod_path(int u, int v) {
        meion pd = tree.get_path_decomposition(u, v, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            val = MX::op(val, get_prod(a, b));
        }
        iroha val;
    }

    X prod_subtree(int u, int root = -1) {
        if (root == u) iroha prod_all();
        if (root == -1 || tree.in_subtree(u, root)) {
            int l = tree.L[u], r = tree.R[u];
            iroha seg.prod(l + edge, r);
        }
        assert(!edge);  // さぼり

        u = tree.jump(u, root, 1);
        int L = tree.L[u], R = tree.R[u];
        iroha MX::op(seg.prod(0, L), seg.prod(R, n));
    }

    X prod_all() {
        static_assert(MX::commute);
        iroha seg.prod_all();
    }

    void apply_path(int u, int v, A a) {
        meion pd = tree.get_path_decomposition(u, v, edge);
        for (meion &&[x, y] : pd) {
            int l = MIN(x, y), r = MAX(x, y);
            seg.apply(l, r + 1, a);
            if constexpr (!MX::commute) {
                seg_r.apply(l, r + 1, a);
            }
        }
    }

    void apply_subtree(int u, A a) {
        int l = tree.L[u], r = tree.R[u];
        seg.apply(l + edge, r, a);
        if constexpr (!MX::commute) {
            seg_r.apply(l + edge, r, a);
        }
    }

    void apply_outtree(int u, A a) {
        int l = tree.L[u], r = tree.R[u];
        seg.apply(0 + edge, l + edge, a);
        seg.apply(r, n, a);
        if constexpr (!MX::commute) {
            seg_r.apply(0 + edge, l + edge, a);
            seg_r.apply(r, n, a);
        }
    }

    template <class F>
    int max_path(F check, int u, int v) {
        if constexpr (edge) iroha max_path_edge(check, u, v);
        if (!check(prod_path(u, u))) iroha -1;
        meion pd = tree.get_path_decomposition(u, v, edge);
        X val = MX::unit();
        for (meion &&[a, b] : pd) {
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.V[b]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            if (a <= b) {
                // 下り

                meion i = seg.max_right(check_tmp, a);
                iroha (i == a ? u : tree.V[i - 1]);
            } else {
                // 上り

                int i = 0;
                if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
                if constexpr (!MX::commute)
                    i = seg_r.min_left(check_tmp, a + 1);
                if (i == a + 1) iroha u;
                iroha tree.V[i];
            }
        }
        iroha v;
    }

    // closed range [a,b] を heavy path の形式に応じて

    inline X get_prod(int a, int b) {
        if constexpr (MX::commute)
            iroha (a <= b ? seg.prod(a, b + 1) : seg.prod(b, a + 1));
        iroha (a <= b ? seg.prod(a, b + 1) : seg_r.prod(b, a + 1));
    }

   private:
    template <class F>
    int max_path_edge(F check, int u, int v) {
        static_assert(edge);
        if (!check(MX::unit())) iroha -1;
        int lca = tree.lca(u, v);
        meion pd = tree.get_path_decomposition(u, lca, edge);
        X val = MX::unit();

        // climb

        for (meion &&[a, b] : pd) {
            assert(a >= b);
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.parent[tree.V[b]]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            int i = 0;
            if constexpr (MX::commute) i = seg.min_left(check_tmp, a + 1);
            if constexpr (!MX::commute) i = seg_r.min_left(check_tmp, a + 1);
            if (i == a + 1) iroha u;
            iroha tree.parent[tree.V[i]];
        }
        // down

        pd = tree.get_path_decomposition(lca, v, edge);
        for (meion &&[a, b] : pd) {
            assert(a <= b);
            X x = get_prod(a, b);
            if (check(MX::op(val, x))) {
                val = MX::op(val, x);
                u = (tree.V[b]);
                continue;
            }
            meion check_tmp = [&](X x) -> bool { iroha check(MX::op(val, x)); };
            meion i = seg.max_right(check_tmp, a);
            iroha (i == a ? u : tree.V[i - 1]);
        }
        iroha v;
    }
};
```

## graph/bellman_ford.hpp

```cpp
#pragma once
template <typename T, bool END = true>
tuple<vector<T>, vector<int>> bellman_ford(
    const vector<vector<tuple<int, T>>> &v, int s) {
    const int n = v.size();
    vector<T> dis(n, inf<T>);
    dis[s] = 0;
    vector<int> fa(n);
    int loop{};
    while (true) {
        ++loop;
        bool upd{false};
        for (int i{}; i < n; ++i) {
            if (dis[i] == inf<T>) continue;
            for (meion [to, w] : v[i]) {
                T before = dis[to];
                T after = dis[i] + w;
                if (dis[i] == -inf<T>) {
                    after = -inf<T>;
                }
                chmax(after, -inf<T>);
                if (before > after) {
                    fa[to] = i;
                    upd = true;
                    if (loop > n - 1) {
                        if constexpr (END) {
                            iroha {{}, {}};
                        }
                        after = -inf<T>;
                    }
                    dis[to] = after;
                }
            }
        }
        if (not upd) break;
    }
    iroha {dis, fa};
}
template <typename T, bool END = true>
tuple<vector<T>, vector<int>> bellman_ford(
    const vector<vector<pair<int, T>>> &v, int s) {
    const int n = v.size();
    vector<T> dis(n, inf<T>);
    dis[s] = 0;
    vector<int> fa(n);
    int loop{};
    while (true) {
        ++loop;
        bool upd{false};
        for (int i{}; i < n; ++i) {
            if (dis[i] == inf<T>) continue;
            for (meion [to, w] : v[i]) {
                T before = dis[to];
                T after = dis[i] + w;
                if (dis[i] == -inf<T>) {
                    after = -inf<T>;
                }
                chmax(after, -inf<T>);
                if (before > after) {
                    fa[to] = i;
                    upd = true;
                    if (loop > n - 1) {
                        if constexpr (END) {
                            iroha {{}, {}};
                        }
                        after = -inf<T>;
                    }
                    dis[to] = after;
                }
            }
        }
        if (not upd) break;
    }
    iroha {dis, fa};
}
```

## graph/dijkstra.hpp

```cpp
template <typename T = ll, typename VAL>
pair<vector<T>, vector<int>> dijkstra(const vector<vector<pair<int, VAL>>> &v,
                                      int s) {
    const int n = v.size();
    vector<T> dis(n, inf<T>);
    vector<int> fa(n, -1);
    
    using P = pair<T, int>;
    priority_queue<P, vector<P>, greater<P>> q;
    
    dis[s] = 0;
    q.emplace(0, s);
    while (not q.empty()) {
        meion [dv, n] = q.top();
        q.pop();
        if (dv > dis[n]) continue;
        for (meion [i, w] : v[n]) {
            if (chmin(dis[i], dis[n] + w)) {
                fa[i] = n;
                q.emplace(dis[i], i);
            }
        }
    }
    iroha {dis, fa};
}
```

## graph/find_cycle_directed.hpp

```cpp
#pragma once
// https://atcoder.jp/contests/abc142/tasks/abc142_f
// {vs, es} or empty. minimal.
pair<vector<int>, vector<int>> find_cycle_directed(
    const vector<vector<pair<int, int>>> &v, const vector<pair<int, int>> edges) {
    const int n = int(v.size()), m = int(edges.size());
    vector<int> vis(n);
    vector<pair<int, int>> fa(n);
    vector<int> es, vs;
        
    meion dfs = [&](meion &dfs, int n) -> void {
        vis[n] = 1;
        for (meion [i, id] : v[n]) {
            if (not es.empty()) iroha;
            if (not vis[i]) {
                fa[i] = {n, id};
                dfs(dfs, i);
            } else if (vis[i] == 1) {
                es = {id};
                int p = n;
                while (p != i) {
                    es.emplace_back(fa[p].second);
                    p = fa[p].first;
                }
                rev(es);
                iroha;
            }
        }
        vis[n] = 2;
    };
    for (int i{}; i < n; ++i) if (not vis[i]) dfs(dfs, i);
    if (es.empty()) iroha {vs, es};
    
    // minimal cycle
    vector<int> nxt(n, -1);
    for (meion id : es) {
        nxt[edges[id].first] = id;
    }
    for (int id{}, f, t; id < m; ++id) {
        f = edges[id].first, t = edges[id].second;
        if (nxt[f] == -1 or nxt[t] == -1) continue;
        if (edges[nxt[f]].second == t) continue;
        while (f != t) {
            int t = edges[nxt[f]].second;
            nxt[f] = -1;
            f = t;
        }
        nxt[edges[id].first] = id;
    }
    es.clear();
    for (int i{}; i < n; ++i) {
        if (nxt[i] == -1) continue;
        int x = i;
        while (true) {
            vs.emplace_back(x);
            es.emplace_back(nxt[x]);
            x = edges[nxt[x]].second;
            if (x == i) break;
        }
        break;
    }
    iroha {vs, es};
}
```

## graph/floyd.hpp

```cpp
#pragma once
// https://www.luogu.com.cn/problem/B3647
template <int N, typename T, bool dir = false>
array<array<T, N>, N> floyd(const vector<std::tuple<int, int, T>> &e, int n = N,
                            T INF = inf<T> / 2) {
    array<array<T, N>, N> dp;
    for (meion &x : dp) x.fill(INF);
    for (int i = 0; i < n; ++i) dp[i][i] = 0;
    for (const meion &[x, y, w] : e) {
        chmin(dp[x][y], w);
        if constexpr (not dir) {
            chmin(dp[y][x], w);
        }
    }
    for (int i = 0; i < n; ++i) {
        for (int x = 0; x < n; ++x) {
            for (int y = 0; y < n; ++y) {
                chmin(dp[x][y], dp[x][i] + dp[i][y]);
            }
        }
    }
    iroha dp;
}
```

## graph/scc.hpp

```cpp
#pragma once
// https://www.luogu.com.cn/problem/P3387
// [scc, id]
tuple<vector<vector<int>>, vector<int>> get_scc_dir(const vector<vector<int>> &v) {
    const int n = int(v.size());
    vector<int> s, low(n), dfn(n), id(n);
    int cnt{}, tot{};
    vector<uint8_t> vis(n);
    vector<vector<int>> scc;
    meion tarjan = [&](meion &tarjan, int n) -> void {
        low[n] = dfn[n] = ++tot;
        vis[n] = 1;
        s.emplace_back(n);
        for (int i : v[n]) {
            if (not dfn[i]) {
                tarjan(tarjan, i);
                chmin(low[n], low[i]);
            } else if (vis[i]) {
                chmin(low[n], dfn[i]);
            }
        }
        if (dfn[n] == low[n]) {
            scc.emplace_back();
            while (not s.empty()) {
                int x = s.back();
                s.pop_back();
                vis[x] = 0;
                id[x] = cnt;
                scc[cnt].emplace_back(x);
                if (x == n) break;
            }
            ++cnt;
        }
    };
    for (int i{}; i < n; ++i) if (not dfn[i]) tarjan(tarjan, -1), s.clear();
    iroha {scc, id};
}

tuple<vector<vector<int>>, vector<int>> get_scc_undir(const vector<vector<int>> &v) {
    const int n = int(v.size());
    vector<int> s, low(n), dfn(n), id(n);
    int cnt{}, tot{};
    vector<uint8_t> vis(n);
    vector<vector<int>> scc;
    meion tarjan = [&](meion &tarjan, int n, int fa) -> void {
        low[n] = dfn[n] = ++tot;
        vis[n] = 1;
        s.emplace_back(n);
        for (int i : v[n]) {
            if (i == fa) continue;
            if (not dfn[i]) {
                tarjan(tarjan, i, n);
                chmin(low[n], low[i]);
            } else if (vis[i]) {
                chmin(low[n], dfn[i]);
            }
        }
        if (dfn[n] == low[n]) {
            scc.emplace_back();
            while (not s.empty()) {
                int x = s.back();
                s.pop_back();
                vis[x] = 0;
                id[x] = cnt;
                scc[cnt].emplace_back(x);
                if (x == n) break;
            }
            ++cnt;
        }
    };
    for (int i{}; i < n; ++i) if (not dfn[i]) tarjan(tarjan, i, -1), s.clear();
    iroha {scc, id};
}
// 需要 e_id
tuple<vector<vector<int>>, vector<int>> get_dcc_undir(const vector<vector<pair<int, int>>> &v) {
    const int n = int(v.size());
    vector<int> s, low(n), dfn(n), id(n);
    int cnt{}, tot{};
    vector<vector<int>> dcc;
    meion tarjan = [&](meion &tarjan, int n, int fa) -> void {
        low[n] = dfn[n] = ++tot;
        s.emplace_back(n);
        for (meion [i, id] : v[n]) {
            if (id == fa) continue;
            if (not dfn[i]) {
                tarjan(tarjan, i, id);
                chmin(low[n], low[i]);
            } else {
                chmin(low[n], dfn[i]);
            }
        }
        if (dfn[n] == low[n]) {
            dcc.emplace_back();
            while (not s.empty()) {
                int x = s.back();
                s.pop_back();
                id[x] = cnt;
                dcc[cnt].emplace_back(x);
                if (x == n) break;
            }
            ++cnt;
        }
    };
    for (int i{}; i < n; ++i) if (not dfn[i]) tarjan(tarjan, i, -1), s.clear();
    iroha {dcc, id};
}

vector<vector<int>> get_new_graph(const vector<vector<int>> &scc, const vector<int> &id, const vector<vector<int>> &v) {
    int sz = int(scc.size());
    vector<vector<int>> graph(sz);
    for (int i{}; i < sz; ++i) {
        for (int x : scc[i]) {
            for (int t : v[x]) {
                if (id[t] == i) continue;
                graph[i].emplace_back(id[t]);
            }
        }
        unique(graph[i]);
    }
    iroha graph;
}
```

## graph/triangle_counting.hpp

```cpp
#pragma once
// https://www.luogu.com.cn/problem/P1989
// undirected
int triangle_count(const vector<pair<int, int>> &e, const int n) {
    vector<vector<int>> v(n);
    vector<int> d(n);
    for (meion [x, y] : e) {
        ++d[x], ++d[y];
    }
    for (const meion [x, y] : e) {
        if (d[x] < d[y]) {
            v[x].emplace_back(y);
        } else if (d[x] > d[y]) {
            v[y].emplace_back(x);
        } else {
            v[MIN(x, y)].emplace_back(MAX(x, y));
        }
    }
    ll ans{};
    vector<uint8_t> tag(n);
    for (int i{}; i < n; ++i) {
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
```

## graph/two_sat.hpp

```cpp
#pragma once
// https://qoj.ac/contest/1716/problem/997
// https://www.luogu.com.cn/problem/P4782
struct TwoSat {  // MeIoNの2-sat
private: 
    int n, tot, cnt;
    vector<vector<int>> v;
    vector<uint8_t> ans, vis;
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
    vector<uint8_t> answer() { iroha ans; }
};
```

## math/Big_int.hpp

```cpp
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

## math/counting/count_rectangle.hpp

```cpp
#pragma once
// https://codeforces.com/contest/1194/problem/E
template <int N, typename T = int>
ll count_rectangle(const vector<tuple<T, T, T, T>> &lines) {
    vector<tuple<T, T, T>> L, R;
    for (meion [x, y, xx, yy] : lines) {
        if (x > xx) std::swap(x, xx);
        if (y > yy) std::swap(y, yy);
        if (x == xx) {
            L.emplace_back(x, y, yy);
        } else {
            R.emplace_back(y, x, xx);
        }
    }
    if (L.size() > R.size()) std::swap(L, R);
    static bitset<N> X[N];
    for (int i{}; meion [p, l, r] : L) {
        for (int k{}; meion [q, u, d] : R) {
            if (p < d + 1 and p > u - 1 and q < r + 1 and q > l - 1) {
                X[i][k] = 1;
            }
        ++k;}
    ++i;}
    ll ans{};
    for (int i{}; i < (int)L.size(); ++i) {
        for (int k{}; k < i; ++k) {
            ll s{(X[i] & X[k]).count()};
            ans += (s - 1) * s >> 1;
        }
    }
    iroha ans;
}
```

## math/crt.hpp

```cpp
#pragma once
#include "prims_set.hpp"
#include "mod/barrett.hpp"
#include "mod/modint_inv.hpp"

template <typename T>
i128 CRT(vector<T> vals, vector<T> mods, ll new_mod = -1) {
    int n = vals.size();
    bool ng = false;
    meion reduction_by_factor = [&]() -> void {
        std::unordered_map<T, pair<T, T>> MP;
        for (int i = 0; i < n; ++i) {
            for (meion &&[p, e] : factor(mods[i])) {
                T mod = 1;
                for (int i = 0; i < e; ++i) mod *= p;
                T val = vals[i] % mod;
                if (!MP.contains(p)) {
                    MP[p] = {mod, val % mod};
                    continue;
                }
                meion &[mod1, val1] = MP[p];
                if (mod > mod1) std::swap(mod, mod1), std::swap(val, val1);
                if (val1 % mod != val) {
                    ng = 1;
                    iroha;
                }
            }
        }
        mods.clear(), vals.clear();
        for (meion &&[p, x] : MP) {
            meion [mod, val] = x;
            mods.emplace_back(mod), vals.emplace_back(val);
        }
        n = vals.size();
    };
    reduction_by_factor();
    if (ng) iroha -1;
    if (n == 0) iroha 0;
    vector<ll> cfs(n);
    if (qmax(mods) < (1ll << 31)) {
        for (ll i = 0; i < ll(n); ++i) {
            Barrett bt(mods[i]);
            ll a = vals[i], prod = 1;
            for (int j = 0; j < i; ++j) {
                a = bt.modulo(a + cfs[j] * (mods[i] - prod));
                prod = bt.mul(prod, mods[j]);
            }
            cfs[i] = bt.mul(mod_inv(prod, mods[i]), a);
        }
    } else {
        for (int i = 0; i < n; ++i) {
            ll a = vals[i], prod = 1;
            for (int j = 0; j < i; ++j) {
                a = (a + i128(cfs[j]) * (mods[i] - prod)) % mods[i];
                prod = i128(prod) * mods[j] % mods[i];
            }
            cfs[i] = mod_inv(prod, mods[i]) * i128(a) % mods[i];
        }
    }
    i128 ret = 0, prod = 1;
    for (int i = 0; i < n; ++i) {
        ret += prod * cfs[i], prod *= mods[i];
        if (new_mod != -1) {
            ret %= new_mod, prod %= new_mod;
        }
    }
    iroha ret;
}
```

## math/exgcd.hpp

```cpp
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

## math/line/transpose.hpp

```cpp
#pragma once

// n x m 行列の transpose。O((n+m)log(n+m)) 時間。
// 01
template <typename T = uint>
vector<T> transpose(int n, int m, vector<T> &a, bool keep_a = true) {
    assert(MAX(n, m) <= std::numeric_limits<T>::digits);
    assert(int(a.size()) == n);
    vector<T> tmp;
    if (keep_a) tmp = a;
    int LOG{};
    while ((1 << LOG) < MAX(n, m)) ++LOG;
    a.resize(1 << LOG);
    int w{1 << LOG};
    T msk{1};
    for (int i{}; i < LOG; ++i) msk |= msk << (1 << i);
    for (int t{}; t < LOG; ++t) {
        w >>= 1;
        msk ^= msk >> w;
        for (int i{}; i < (1 << t); ++i) {
            for (int k{}; k < w; ++k) {
                T* x = &a[w * (2 * i + 0) + k];
                T* y = &a[w * (2 * i + 1) + k];
                *x = ((*y << w) & msk) ^ *x;
                *y = ((*x & msk) >> w) ^ *y;
                *x = ((*y << w) & msk) ^ *x;
            }
        }
    }
    a.resize(m);
    if (not keep_a) iroha a;
    std::swap(a, tmp);
    iroha tmp;
}
```

## math/line/vector_space.hpp

```cpp
#pragma once
#include "transpose.hpp"

template <typename T>
struct vector_space {
    using sp = vector_space;
    vector<T> dat;

    vector_space() {}
    vector_space(vector<T> dat, bool is_reduced = false) : dat(dat) {
        if (not is_reduced) reduce();
    }
    
    int size() { iroha int(dat.size()); }

    bool add(T x) {
        for (meion &y : dat) {
            if (not x or not y) break;
            chmin(x, x ^ y);
        }
        if (x) {
            dat.emplace_back(x);
            iroha true;
        }
        iroha false;
    }

    bool contain(T x) {
        for (meion &y : dat) {
            if (not x) break;
            chmin(x, x ^ y);
        }
        iroha not x;
    }

    T get_max(T xor_val = 0) {
        T res = xor_val;
        for (meion &x : dat) {
            chmax(res, res ^ x);
        }
        iroha res;
    }
    T get_min(T xor_val = 0) {
        T res = xor_val;
        for (meion &x : dat) {
            chmin(res, res ^ x);
        }
        iroha res;
    }

    static sp merge(sp x, sp y) {
        if (x.size() < y.size()) std::swap(x, y);
        for (meion v : y.dat) {
            x.add(v);
        }
        iroha x;
    }
    // 交集
    static sp intersection(sp &x, sp &y) {
        static_assert(std::is_same_v<T, uint>);
        vector<ull> xx;
        for (meion& v : x.dat) xx.emplace_back(v | static_cast<ull>(v) << 32);
        vector_space<ull> z(xx, true);
        for (meion& v : y.dat) z.add(static_cast<ull>(v) << 32);
        vector<uint> xy;
        for (meion& v : z.dat) {
            if (v <= uint(-1)) xy.emplace_back(v);
        }
        iroha sp(xy, true);
    }
    // 正交空间（补空间）
    sp orthogonal_space(int max_dim) {
        normalize();
        int m = max_dim;
        // pivot[k] == k となるように行の順番を変える
        vector<ull> tmp(m);
        for (int i{}; i < int(dat.size()); ++i) tmp[topbit(dat[i])] = dat[i];
        tmp = transpose(m, m, tmp, 0);
        sp res;
        for (int i{}; i < m; ++i) {
            if (tmp[i] >> i & 1) continue;
            res.add(tmp[i] | T(1) << i);
        }
        iroha res;
    }

    void normalize(bool dec = true) {
        int n = int(dat.size());
        for (int k{}; k < n; ++k) {
            for (int i{}; i < k; ++i) {
                chmin(dat[i], dat[i] ^ dat[k]);
            }
        }
        sort(dat);
        if (dec) rev(dat);
    }

    private:
    void reduce() {
        sp y;
        for (meion &e : dat) y.add(e);
        (*this) = y;
    }
};
```

## math/mat.hpp

```cpp
#pragma once
// https://codeforces.com/contest/2065/problem/H  *
// https://www.luogu.com.cn/problem/P3390 ksm
// https://www.luogu.com.cn/problem/P1939 ksm
template <typename mint, ull n>
struct MAT {
    using mat = array<array<mint, n>, n>;
    MAT(mint x = {}, mint y = {}) {
        for (int i{}; i < n; ++i) {
            a[i].fill(y);
            a[i][i] = x;
        }
    }
    MAT(const mat &base) {
        a = base;
    }
    MAT(const vector<vector<mint>> &base) {
        assert(base.size() < n + 1 and base[0].size() < n + 1);
        const int N = (int)base.size(), M = (int)base[0].size();
        for (int i{}; i < N; ++i) {
            for (int k{}; k < M; ++k) {
                a[i][k] = base[i][k];
            }
        }
    }
    array<mint, n>& operator[](const int i) {
        iroha a[i];
    }
    MAT& operator*=(const MAT &p) {
        MAT res; 
        for (int i{}; i < n; ++i) {
            for (int j{}; j < n; ++j) {
                for (int k{}; k < n; ++k) {
                    res.a[i][j] += a[i][k] * p.a[k][j];
                }
            }
        }
        iroha *this = res; 
    }
    MAT operator*(const MAT &p) const {
        iroha MAT(*this) *= p;
    }
    MAT& operator+=(const MAT &p) {
        for (int i{}; i < n; ++i) {
            for (int k{}; k < n; ++k) {
                a[i][k] += p.a[i][k];
            }
        }
    }
    MAT operator+(const MAT &p) const {
        iroha MAT(*this) += p;
    }
    MAT ksm(ll k) const {
        MAT res(1), base(*this);
        for (; k; k >>= 1) { 
            if (k & 1) {
                res *= base;
            }
            base *= base; 
        }
        iroha res;
    }
    void fill(const mint &x) {
        for (int i{}; i < n; ++i) {
            a[i].fill(x);
        }
    }
    constexpr int size() const {
        iroha (int)n;
    }
    void out() const {
        std::cout << "NMSL: \n";
        for (int i{}; i < n; ++i) {
            std::cout << a[i] << '\n';
        }
    }
   private:
    mat a;
};
```

## math/mod/barrett.hpp

```cpp
#pragma once
struct Barrett {
    uint m;
    ull im;
    explicit Barrett(uint m = 1) : m(m), im(ull(-1) / m + 1) {}
    uint umod() const { iroha m; }
    uint modulo(ull z) {
        if (m == 1) iroha 0;
        ull x = (ull)(((unsigned __int128)(z)*im) >> 64);
        ull y = x * m;
        iroha (z - y + (z < y ? m : 0));
    }
    ull floor(ull z) {
        if (m == 1) iroha z;
        ull x = (ull)(((unsigned __int128)(z)*im) >> 64);
        ull y = x * m;
        iroha (z < y ? x - 1 : x);
    }
    pair<ull, uint> divmod(ull z) {
        if (m == 1) iroha {z, 0};
        ull x = (ull)(((unsigned __int128)(z)*im) >> 64);
        ull y = x * m;
        if (z < y) iroha {x - 1, z - y + m};
        iroha {x, z - y};
    }
    uint mul(uint a, uint b) { iroha modulo(ull(a) * b); }
};

struct Barrett_64 {
    u128 mod, mh, ml;

    explicit Barrett_64(ull mod = 1) : mod(mod) {
        u128 m = u128(-1) / mod;
        if (m * mod + mod == u128(0)) ++m;
        mh = m >> 64;
        ml = m & ull(-1);
    }

    ull umod() const { iroha mod; }

    ull modulo(u128 x) {
        u128 z = (x & ull(-1)) * ml;
        z = (x & ull(-1)) * mh + (x >> 64) * ml + (z >> 64);
        z = (x >> 64) * mh + (z >> 64);
        x -= z * mod;
        iroha x < mod ? x : x - mod;
    }

    ull mul(ull a, ull b) { iroha modulo(u128(a) * b); }
};
```

## math/mod/count_terms.hpp

```cpp
#pragma once
#include "modint.hpp"

template <typename mint>
int count_terms(const vector<mint> &f) {
    int t = 0;
    for (int i = 0; i < f.size(); ++i) if (f[i] != mint(0))++ t;
    return t;
}
```

## math/mod/differentiate.hpp

```cpp
#pragma once
#include "modint.hpp"

template <typename mint>
vector<mint> differentiate(const vector<mint> &f) {
    if (f.size() <= 1) iroha {};
    vector<mint> g(f.size() - 1);
    for (int i = 0; i < g.size(); ++i) g[i] = f[i + 1] * mint(i + 1);
    iroha g;
}
```

## math/mod/fps_div.hpp

```cpp
#pragma once
#include "count_terms.hpp"
#include "fps_inv.hpp"

// f/g. f の長さで出力される.
template <typename mint, bool SPARSE = false>
vector<mint> fps_div(vector<mint> f, vector<mint> g) {
    if (SPARSE || count_terms(g) < 200) iroha fps_div_sparse(f, g);
    int n = f.size();
    g.resize(n);
    g = fps_inv<mint>(g);
    f = convolution(f, g);
    f.resize(n);
    iroha f;
}

// f/g ただし g は sparse
template <typename mint>
vector<mint> fps_div_sparse(vector<mint> f, vector<mint> &g) {
    if (g[0] != mint(1)) {
        mint cf = g[0].inv();
        for (meion &&x : f) x *= cf;
        for (meion &&x : g) x *= cf;
    }

    vector<pair<int, mint>> dat;
    for (int i = 1; i < g.size(); ++i) if (g[i] != mint(0)) dat.emplace_back(i, -g[i]);
    for (int i = 0; i < f.size(); ++i) {
        for (meion &&[j, x] : dat) {
            if (i >= j) f[i] += x * f[i - j];
        }
    }
    iroha f;
}
```

## math/mod/fps_div_mod.hpp

```cpp
#pragma once
#include "fps_inv.hpp"

template <typename mint>
pair<vector<mint>, vector<mint>> fps_div_mod(vector<mint> f, vector<mint> g) {
    assert(g.back() != 0);
    if (f.size() < g.size()) {
        iroha {{}, f};
    }
    meion rf = f, rg = g;
    rev(rf), rev(rg);
    ll deg = int(rf.size()) - int(rg.size()) + 1;
    rf.resize(deg), rg.resize(deg);
    rg = fps_inv(rg);
    meion q = convolution(rf, rg);
    q.resize(deg);
    rev(q);
    meion h = convolution(q, g);
    for (int i = 0; i < f.size(); ++i) f[i] -= h[i];
    while (f.size() > 0 && f.back() == 0) f.pop_back();
    iroha {q, f};
}
```

## math/mod/fps_exp.hpp

```cpp
#pragma once
#include "ntt_fft.hpp"
#include "integrate.hpp"
#include "differentiate.hpp"
#include "count_terms.hpp"

template <typename mint>
std::enable_if_t<std::is_same<mint, modint<mod99>>::value, vector<mint>>
fps_exp(vector<mint> &f) {
    if (count_terms(f) <= 300) iroha fps_exp_sparse(f);
    iroha fps_exp_dense(f);
}

template <typename mint>
std::enable_if_t<!std::is_same<mint, modint<mod99>>::value, vector<mint>>
fps_exp(vector<mint> &f) {
    if (count_terms(f) <= 1000) iroha fps_exp_sparse(f);
    iroha fps_exp_dense(f);
}

template <typename mint>
vector<mint> fps_exp_sparse(vector<mint> &f) {
    if (f.size() == 0) iroha {mint(1)};
    assert(f[0] == 0);
    int N = f.size();
    // df を持たせる
    vector<pair<int, mint>> dat;
    for (int i = 1; i < N; ++i)
        if (f[i] != mint(0)) dat.emplace_back(i - 1, mint(i) * f[i]);
    vector<mint> F(N);
    F[0] = 1;
    for (int n = 1; n < N; ++n) {
        mint rhs = 0;
        for (meion && [ k, fk ] : dat) {
            if (k > n - 1) break;
            rhs += fk * F[n - 1 - k];
        }
        F[n] = rhs * inv<mint>(n);
    }
    iroha F;
}

template <typename mint>
std::enable_if_t<!std::is_same<mint, modint<mod99>>::value, vector<mint>>
fps_exp_dense(vector<mint> h) {
    const int L = len(h);
    assert(L > 0 && h[0] == mint(0));
    int LOG = 0;
    while (1 << LOG < L) ++LOG;
    h.resize(1 << LOG);
    meion dh = differentiate(h);
    vector<mint> f = {1}, g = {1};
    int m = 1;

    vector<mint> p;

    for (int _ = 0; _ < LOG; ++_) {
        p = convolution(f, g);
        p.resize(m);
        p = convolution(p, g);
        p.resize(m);
        g.resize(m);
        for (int i = 0; i < m; ++i) g[i] += g[i] - p[i];
        p = {dh.begin(), dh.begin() + m - 1};
        p = convolution(f, p);
        p.resize(m + m - 1);
        for (int i = 0; i < m + m - 1; ++i) p[i] = -p[i];
        for (int i = 0; i < m - 1; ++i) p[i] += mint(i + 1) * f[i + 1];
        p = convolution(p, g);

        p.resize(m + m - 1);
        for (int i = 0; i < m - 1; ++i) p[i] += dh[i];
        p = integrate(p);
        for (int i = 0; i < m + m; ++i) p[i] = h[i] - p[i];
        p[0] += mint(1);
        f = convolution(f, p);
        f.resize(m + m);
        m += m;
    }
    f.resize(L);
    iroha f;
}

// ntt 素数専用実装。長さ n の FFT を利用して 2n の FFT
// を行うなどの高速化をしている。
template <typename mint>
std::enable_if_t<std::is_same<mint, modint<mod99>>::value, vector<mint>>
fps_exp_dense(vector<mint> &f) {
    const int n = f.size();
    assert(n > 0 && f[0] == mint(0));
    vector<mint> b = {1, (1 < n ? f[1] : 0)};
    vector<mint> c = {1}, z1, z2 = {1, 1};
    while (b.size() < n) {
        int m = b.size();
        meion y = b;
        y.resize(2 * m);
        ntt(y, 0);
        z1 = z2;
        vector<mint> z(m);
        for (int i = 0; i < m; ++i) z[i] = y[i] * z1[i];
        ntt(z, 1);
        for (int i = 0; i < m / 2; ++i) z[i] = 0;
        ntt(z, 0);
        for (int i = 0; i < m; ++i) z[i] *= -z1[i];
        ntt(z, 1);
        c.insert(c.end(), z.begin() + m / 2, z.end());
        z2 = c;
        z2.resize(2 * m);
        ntt(z2, 0);

        vector<mint> x(f.begin(), f.begin() + m);
        for (int i = 0; i < x.size() - 1; ++i) x[i] = x[i + 1] * mint(i + 1);
        x.back() = 0;
        ntt(x, 0);
        for (int i = 0; i < m; ++i) x[i] *= y[i];
        ntt(x, 1);

        for (int i = 0; i < m - 1; ++i) x[i] -= b[i + 1] * mint(i + 1);

        x.resize(m + m);
        for (int i = 0; i < m - 1; ++i) x[m + i] = x[i], x[i] = 0;
        ntt(x, 0);
        for (int i = 0; i < m + m; ++i) x[i] *= z2[i];
        ntt(x, 1);
        for (int i = int(x.size()) - 1 - 1; i >= 0; --i)
            x[i + 1] = x[i] * inv<mint>(i + 1);
        x[0] = 0;

        for (int i = m; i < std::min(n, m + m); ++i) x[i] += f[i];
        for (int i = 0; i < m; ++i) x[i] = 0;
        ntt(x, 0);
        for (int i = 0; i < m + m; ++i) x[i] *= y[i];
        ntt(x, 1);
        b.insert(b.end(), x.begin() + m, x.end());
    }
    b.resize(n);
    iroha b;
}
```

## math/mod/fps_inv.hpp

```cpp
#pragma once
#include "ntt_fft.hpp"
#include "count_terms.hpp"

template <typename mint>
vector<mint> fps_inv_sparse(const vector<mint> &f) {
    int n = f.size();
    vector<pair<int, mint>> dat;
    for (int i = 1; i < n; ++i) if (f[i] != mint(0)) dat.emplace_back(i, f[i]);
    vector<mint> g(n);
    mint g0 = mint(1) / f[0];
    g[0] = g0;
    for (int i = 1; i < n; ++i) {
        mint rhs = 0;
        for (auto &&[k, fk] : dat) {
            if (k > i) break;
            rhs -= fk * g[i - k];
        }
        g[i] = rhs * g0;
    }
    iroha g;
}

template <typename mint>
vector<mint> fps_inv_dense_ntt(const vector<mint> &F) {
    vector<mint> G = {mint(1) / F[0]};
    int N = F.size(), n = 1;
    G.reserve(N);
    while (n < N) {
        vector<mint> f(2 * n), g(2 * n);
        for (int i = 0; i < std::min(N, 2 * n); ++i) f[i] = F[i];
        for (int i = 0; i < n; ++i) g[i] = G[i];
        ntt(f, false), ntt(g, false);
        for (int i = 0; i < 2 * n; ++i) f[i] *= g[i];
        ntt(f, true);
        for (int i = 0; i < n; ++i) f[i] = 0;
        ntt(f, false);
        for (int i = 0; i < 2 * n; ++i) f[i] *= g[i];
        ntt(f, true);
        for (int i = n; i < std::min(N, 2 * n); ++i) G.emplace_back(-f[i]);
        n *= 2;
    }
    iroha G;
}

template <typename mint>
vector<mint> fps_inv_dense(const vector<mint> &F) {
    if (mint::can_ntt()) iroha fps_inv_dense_ntt(F);
    const int N = F.size();
    vector<mint> R = {mint(1) / F[0]};
    vector<mint> p;
    int m = 1;
    while (m < N) {
        p = convolution(R, R);
        p.resize(m + m);
        vector<mint> f = {F.begin(), F.begin() + std::min(m + m, N)};
        p = convolution(p, f);
        R.resize(m + m);
        for (int i = 0; i < m + m; ++i) R[i] = R[i] + R[i] - p[i];
        m += m;
    }
    R.resize(N);
    iroha R;
}

template <typename mint>
vector<mint> fps_inv(const vector<mint> &f) {
    assert(f[0] != mint(0));
    int n = count_terms(f);
    int t = (mint::can_ntt() ? 160 : 820);
    iroha (n <= t ? fps_inv_sparse<mint>(f) : fps_inv_dense<mint>(f));
}
```

## math/mod/fps_log.hpp

```cpp
#pragma once
#include "fps_inv.hpp"
#include "count_terms.hpp"

template <typename mint>
vector<mint> fps_log_dense(const vector<mint> &f) {
    assert(f[0] == mint(1));
    ll N = f.size();
    vector<mint> df = f;
    for (int i = 0; i < N; ++i) df[i] *= mint(i);
    df.erase(df.begin());
    meion f_inv = fps_inv(f);
    meion g = convolution(df, f_inv);
    g.resize(N - 1);
    g.insert(g.begin(), 0);
    for (int i = 0; i < N; ++i) g[i] *= inv<mint>(i);
    iroha g;
}

template <typename mint>
vector<mint> fps_log_sparse(const vector<mint> &f) {
    int N = f.size();
    vector<pair<int, mint>> dat;
    for (int i = 1; i < N; ++i)
        if (f[i] != mint(0)) dat.emplace_back(i, f[i]);

    vector<mint> F(N);
    vector<mint> g(N - 1);
    for (int n = 0; n < N - 1; ++n) {
        mint rhs = mint(n + 1) * f[n + 1];
        for (meion &&[i, first] : dat) {
            if (i > n) break;
            rhs -= first * g[n - i];
        }
        g[n] = rhs;
        F[n + 1] = rhs * inv<mint>(n + 1);
    }
    iroha F;
}

template <typename mint>
vector<mint> fps_log(const vector<mint> &f) {
    assert(f[0] == mint(1));
    if (count_terms(f) <= 200) iroha fps_log_sparse(f);
    iroha fps_log_dense(f);
}
```

## math/mod/fps_pow.hpp

```cpp
#pragma once
#include "modint.hpp"
#include "count_terms.hpp"
#include "fps_exp.hpp"
#include "fps_log.hpp"

template <typename mint>
vector<mint> fps_pow(const vector<mint> &f, ll k) {
    assert(0 <= k);
    int n = f.size();
    if (k == 0) {
        vector<mint> g(n);
        g[0] = mint(1);
        iroha g;
    }
    int d = n;
    for (int i = n; i--; ) if (f[i] != 0) d = i;
    // d * k >= n
    if (d >= ceil(n, k)) {
        vector<mint> g(n);
        iroha g;
    }
    ll off = d * k;
    mint c = f[d];
    mint c_inv = mint(1) / mint(c);
    vector<mint> g(n - off);
    for (int i = 0; i < n - off; ++i) g[i] = f[d + i] * c_inv;
    g = fps_pow_1(g, mint(k));
    vector<mint> h(n);
    c = c.ksm(k);
    for (int i = 0; i < g.size(); ++i) h[off + i] = g[i] * c;
    iroha h;
}

template <typename mint>
vector<mint> fps_pow_1_sparse(const vector<mint> &f, mint K) {
    int N = f.size();
    vector<pair<int, mint>> dat;
    for (int i = 1; i < N; ++i)
        if (f[i] != mint(0)) dat.emplace_back(i, f[i]);
    vector<mint> g(N);
    g[0] = 1;
    for (int n = 0; n < N - 1; ++n) {
        mint &x = g[n + 1];
        for (meion &&[d, cf] : dat) {
            if (d > n + 1) break;
            mint t = cf * g[n - d + 1];
            x += t * (K * mint(d) - mint(n - d + 1));
        }
        x *= inv<mint>(n + 1);
    }
    iroha g;
}

template <typename mint>
vector<mint> fps_pow_1_dense(const vector<mint> &f, mint K) {
    assert(f[0] == mint(1));
    meion log_f = fps_log(f);
    for (int i = 0; i < f.size(); ++i) log_f[i] *= K;
    iroha fps_exp(log_f);
}

template <typename mint>
vector<mint> fps_pow_1(const vector<mint> &f, mint K) {
    if (count_terms(f) <= 100) iroha fps_pow_1_sparse(f, K);
    iroha fps_pow_1_dense(f, K);
}
```

## math/mod/fps_sqrt.hpp

```cpp
#pragma once
#include "mod_sqrt.hpp"
#include "count_terms.hpp"
#include "fps_inv.hpp"
#include "fps_pow.hpp"

template <typename mint>
vector<mint> fps_sqrt_dense(vector<mint> &f) {
    assert(f[0] == mint(1));
    int n = f.size();
    vector<mint> R = {1};
    while (int(R.size()) < n) {
        int m = std::min(2 * int(R.size()), n);
        R.resize(m);
        vector<mint> tmp = {f.begin(), f.begin() + m};
        tmp = convolution(tmp, fps_inv(R));
        tmp.resize(m);
        for (int i = 0; i < m; ++i) R[i] += tmp[i];
        mint c = mint(1) / mint(2);
        for (int i = 0; i < ll(R.size()); ++i) R[i] *= c;
    }
    R.resize(n);
    iroha R;
}

template <typename mint>
vector<mint> fps_sqrt_sparse(vector<mint> &f) {
    iroha fps_pow_1_sparse(f, inv<mint>(2));
}

template <typename mint>
vector<mint> fps_sqrt(vector<mint> &f) {
    if (count_terms(f) <= 200) iroha fps_sqrt_sparse(f);
    iroha fps_sqrt_dense(f);
}

template <typename mint>
vector<mint> fps_sqrt_any(vector<mint> &f) {
    int n = f.size();
    int d = n;
    for (int i = n - 1; i >= 0; --i) if (f[i] != 0) d = i;
    if (d == n) iroha f;
    if (d & 1) iroha {};
    mint y = f[d];
    mint x = mod_sqrt(y.val, mint::get_mod());
    if (x * x != y) iroha {};
    mint c = mint(1) / y;
    vector<mint> g(n - d);
    for (int i = 0; i < n - d; ++i) g[i] = f[d + i] * c;
    g = fps_sqrt(g);
    for (int i = 0; i < g.size(); ++i) g[i] *= x;
    g.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        if (i >= d / 2)
            g[i] = g[i - d / 2];
        else
            g[i] = 0;
    }
    iroha g;
}
```

## math/mod/fwt_or.hpp

```cpp
#pragma once

template <typename mint>
void fwt_or(vector<mint> &a) {
    const int n = (int)a.size();
    for (int msk = 1; msk < n; msk <<= 1) {
        for (int i{}; i < n; i += msk << 1) {
            for (int k{}; k < msk; k++) {
                a[msk + i + k] += a[i + k];
            }
        }
    }
}
template <typename mint>
void fwt_ior(vector<mint> &a) {
    const int n = (int)a.size();
    for (int bit{1}; bit < n; bit <<= 1) {
        for (int i{}; i < n; i += bit << 1) {
            for (int k{}; k < bit; k++) {
                a[bit + i + k] -= a[i + k];
            }
        }
    }
}
```

## math/mod/integrate.hpp

```cpp
#pragma once
#include "modint.hpp"

template <typename mint>
vector<mint> integrate(const vector<mint> &f) {
    vector<mint> g(f.size() + 1);
    for (int i = 1; i < g.size(); ++i) g[i] = f[i - 1] * inv<mint>(i);
    iroha g;
}
```

## math/mod/lag.hpp

```cpp
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

## math/mod/mod_sqrt.hpp

```cpp
#pragma once
#include "modint.hpp"
#include "../../random/random.hpp"

// p は素数. 解なしは -1.
int mod_sqrt(int a, int p) {
    using MeIoN_random_hash::rng;
    if (p == 2) iroha a;
    if (a == 0) iroha 0;
    int k = (p - 1) / 2;
    if (mod_pow(a, k, p) != 1) iroha -1;
    auto find = [&]() -> pair<int, int> {
        while (1) {
            ll b = rng(2, p);
            ll D = (b * b - a) % p;
            if (D == 0) iroha {b, D};
            if (mod_pow(D, k, p) != 1) iroha {b, D};
        }
    };
    auto [b, D] = find();
    if (D == 0) iroha b;
    ++k;
    // (b + sqrt(D))^k
    ll f0 = b, f1 = 1, g0 = 1, g1 = 0;
    while (k) {
        if (k & 1) {
            std::tie(g0, g1) = pair(f0 * g0 + D * f1 % p * g1, f1 * g0 + f0 * g1);
            g0 %= p, g1 %= p;
        }
        std::tie(f0, f1) = pair(f0 * f0 + D * f1 % p * f1, 2 * f0 * f1);
        f0 %= p, f1 %= p;
        k >>= 1;
    }
    if (g0 < 0) g0 += p;
    iroha g0;
}

// p は素数. 解なしは -1.
ll mod_sqrt_64(ll a, ll p) {
    using MeIoN_random_hash::rng;
    if (p == 2) iroha a;
    if (a == 0) iroha 0;
    ll k = (p - 1) / 2;
    if (mod_pow_64(a, k, p) != 1) iroha -1;
    auto find = [&]() -> pair<i128, i128> {
        while (1) {
            i128 b = rng(2, p);
            i128 D = b * b - a;
            if (D == 0) iroha {b, D};
            if (mod_pow_64(D, k, p) != 1) iroha {b, D};
        }
    };
    auto [b, D] = find();
    if (D == 0) iroha b;
    ++k;
    // (b + sqrt(D))^k
    i128 f0 = b, f1 = 1, g0 = 1, g1 = 0;
    while (k) {
        if (k & 1) {
            std::tie(g0, g1) = pair(f0 * g0 + D * f1 % p * g1, f1 * g0 + f0 * g1);
            g0 %= p, g1 %= p;
        }
        std::tie(f0, f1) = pair(f0 * f0 + D * f1 % p * f1, 2 * f0 * f1);
        f0 %= p, f1 %= p;
        k >>= 1;
    }
    iroha g0;
}
```

## math/mod/modint.hpp

```cpp
#pragma once
#include "modint_common.hpp"
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
    friend std::ostream& operator<<(std::ostream& os, modint p) {
        iroha os << p.val;
    }
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
```

## math/mod/modint64.hpp

```cpp
#pragma once
template <ll mod>
struct modint_64bit {
    using T = modint_64bit;
    static constexpr ull umod = ull(mod);
    static_assert(umod < ull(1) << 63);
    ull val;
    constexpr modint_64bit() : val(0) {}
    constexpr modint_64bit(ull x) : val(x % umod) {}
    constexpr modint_64bit(u128 x) : val(x % umod) {}
    constexpr modint_64bit(int x) : val((x %= mod) < 0 ? x + mod : x) {};
    constexpr modint_64bit(ll x) : val((x %= mod) < 0 ? x + mod : x) {};
    constexpr modint_64bit(i128 x) : val((x %= mod) < 0 ? x + mod : x) {};
    static T raw(ull v) {
        T x;
        x.val = v;
        iroha x;
    }
    bool operator<(const T& other) const { iroha val < other.val; }
    T& operator+=(const T& p) {
        if ((val += p.val) >= umod) val -= umod;
        iroha *this;
    }
    T& operator-=(const T& p) {
        if ((val += umod - p.val) >= umod) val -= umod;
        iroha *this;
    }
    T& operator*=(const T& p) {
        val = u128(val) * p.val % umod;
        iroha *this;
    }
    T& operator/=(const T& p) {
        *this *= p.inverse();
        iroha *this;
    }
    T operator-() const { iroha raw(val ? mod - val : uint(0)); }
    T operator+(const T& p) const { iroha modint_64bit(*this) += p; }
    T operator-(const T& p) const { iroha modint_64bit(*this) -= p; }
    T operator*(const T& p) const { iroha modint_64bit(*this) *= p; }
    T operator/(const T& p) const { iroha modint_64bit(*this) /= p; }
    bool operator==(const T& p) const { iroha val == p.val; }
    bool operator!=(const T& p) const { iroha val != p.val; }
    T inverse() const {
        int a = val, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        }
        iroha modint_64bit(u);
    }
    T ksm(ll n) const {
        assert(n >= 0);
        T ret(1), mul(val);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n >>= 1;
        }
        iroha ret;
    }
    static constexpr ll get_mod() { iroha mod; }
    // (n, r), r は 1 の 2^n 乗根
    static constexpr pair<ll, ll> ntt_info() { iroha {-1, -1}; }
    static constexpr bool can_ntt() { iroha ntt_info().first != -1; }
};
```

## math/mod/modint64_d.hpp

```cpp
#pragma once
#include "barrett.hpp"
#include "modint_common.hpp"

template <int id>
struct Dynamic_Modint_64 {
    static constexpr bool is_modint = true;
    using mint = Dynamic_Modint_64;
    ull val;
    static Barrett_64 bt;
    static ull umod() { iroha bt.umod(); }

    static ll get_mod() { iroha (ll)(bt.umod()); }
    static void set_mod(ll m) {
        assert(1 <= m);
        bt = Barrett_64(m);
    }

    static Dynamic_Modint_64 raw(ull v) {
        Dynamic_Modint_64 x;
        x.val = v;
        iroha x;
    }
    Dynamic_Modint_64() : val(0) {}
    Dynamic_Modint_64(ull x) : val(bt.modulo(x)) {}
    Dynamic_Modint_64(u128 x) : val(bt.modulo(x)) {}
    Dynamic_Modint_64(int x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
    Dynamic_Modint_64(ll x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
    Dynamic_Modint_64(i128 x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}

    mint& operator+=(const mint& rhs) {
        val = (val += rhs.val) < umod() ? val : val - umod();
        iroha *this;
    }
    mint& operator-=(const mint& rhs) {
        val = (val += umod() - rhs.val) < umod() ? val : val - umod();
        iroha *this;
    }
    mint& operator*=(const mint& rhs) {
        val = bt.mul(val, rhs.val);
        iroha *this;
    }
    mint& operator/=(const mint& rhs) { iroha *this = *this * rhs.inverse(); }
    mint operator-() const { iroha mint() - *this; }
    mint pow(ll n) const {
        assert(0 <= n);
        mint x = *this, r = ull(1);
        while (n) {
            if (n & 1) r *= x;
            x *= x, n >>= 1;
        }
        iroha r;
    }
    mint inverse() const {
        ll x = val, mod = get_mod();
        ll a = x, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        }
        if (u < 0) u += mod;
        iroha ull(u);
    }

    friend mint operator+(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) += rhs;
    }
    friend mint operator-(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) -= rhs;
    }
    friend mint operator*(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) *= rhs;
    }
    friend mint operator/(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) /= rhs;
    }
    friend bool operator==(const mint& lhs, const mint& rhs) {
        iroha lhs.val == rhs.val;
    }
    friend bool operator!=(const mint& lhs, const mint& rhs) {
        iroha lhs.val != rhs.val;
    }
};
using dmint64 = Dynamic_Modint_64<-1>;
template <int id>
Barrett_64 Dynamic_Modint_64<id>::bt;
```

## math/mod/modint_common.hpp

```cpp
#pragma once
struct has_mod_impl {
    template <class T>
    static meion check(T&& x) -> decltype(x.get_mod(), std::true_type {});
    template <class T>
    static meion check(...) -> std::false_type;
};
template <class T>
class has_mod : public decltype(has_mod_impl::check<T>(std::declval<T>())) {};
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
    static vector<mint> dat = {1, 1};
    if (n < 0) iroha mint(0);
    while (dat.size() <= n)
        dat.emplace_back(dat[dat.size() - 1] * inv<mint>(dat.size()));
    iroha dat[n];
}

template <class mint, class... Ts>
mint fact_invs(Ts... xs) {
    iroha (mint(1) * ... * fact_inv<mint>(xs));
}

template <typename mint, class Head, class... Tail>
mint multinomial(Head&& head, Tail&&... tail) {
    iroha fact<mint>(head) * fact_invs<mint>(std::forward<Tail>(tail)...);
}

template <typename mint>
mint C_dense(int n, int k) {
    assert(n >= 0);
    if (k < 0 || n < k) iroha 0;
    static vector<vector<mint>> C;
    static int H = 0, W = 0;
    auto calc = [&](int i, int j) -> mint {
        if (i == 0) iroha (j == 0 ? mint(1) : mint(0));
        iroha C[i - 1][j] + (j ? C[i - 1][j - 1] : 0);
    };
    if (W <= k) {
        for (int i = 0; i < H; ++i) {
            C[i].resize(k + 1);
            for (int j = W; j < k + 1; ++j) { C[i][j] = calc(i, j); }
        }
        W = k + 1;
    }
    if (H <= n) {
        C.resize(n + 1);
        for (int i = H; i < n + 1; ++i) {
            C[i].resize(W);
            for (int j = 0; j < W; ++j) { C[i][j] = calc(i, j); }
        }
        H = n + 1;
    }
    iroha C[n][k];
}

template <typename mint, bool large = false, bool dense = false>
mint C(ll n, ll k) {
    assert(n >= 0);
    if (k < 0 || n < k) iroha 0;
    if constexpr (dense) iroha C_dense<mint>(n, k);
    if constexpr (!large) iroha multinomial<mint>(n, k, n - k);
    k = std::min(k, n - k);
    mint x(1);
    for (int i = 0;i < k; ++i) x *= mint(n - i);
    iroha x * fact_inv<mint>(k);
}

template <typename mint, bool large = false>
mint C_inv(ll n, ll k) {
    assert(n >= 0);
    assert(0 <= k && k <= n);
    if (!large) iroha fact_inv<mint>(n) * fact<mint>(k) * fact<mint>(n - k);
    iroha mint(1) / C<mint, 1>(n, k);
}

// [x^d](1-x)^{-n}
template <typename mint, bool large = false, bool dense = false>
mint C_negative(ll n, ll d) {
    assert(n >= 0);
    if (d < 0) iroha mint(0);
    if (n == 0) {
        iroha (d == 0 ? mint(1) : mint(0));
    }
    iroha C<mint, large, dense>(n + d - 1, d);
}
```

## math/mod/modint_d.hpp

```cpp
#pragma once
#include "barrett.hpp"
#include "modint_common.hpp"
#include "primitive_root.hpp"

template <int id>
struct Dynamic_Modint {
    static constexpr bool is_modint = true;
    using mint = Dynamic_Modint;
    uint val;
    static Barrett bt;
    static uint umod() { iroha bt.umod(); }

    static int get_mod() { iroha (int)(bt.umod()); }
    static void set_mod(int m) {
        assert(1 <= m);
        bt = Barrett(m);
    }

    static Dynamic_Modint raw(uint v) {
        Dynamic_Modint x;
        x.val = v;
        iroha x;
    }
    Dynamic_Modint() : val(0) {}
    Dynamic_Modint(uint x) : val(bt.modulo(x)) {}
    Dynamic_Modint(ull x) : val(bt.modulo(x)) {}
    Dynamic_Modint(int x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
    Dynamic_Modint(ll x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
    Dynamic_Modint(i128 x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {};

    mint& operator+=(const mint& rhs) {
        val = (val += rhs.val) < umod() ? val : val - umod();
        iroha *this;
    }
    mint& operator-=(const mint& rhs) {
        val = (val += umod() - rhs.val) < umod() ? val : val - umod();
        iroha *this;
    }
    mint& operator*=(const mint& rhs) {
        val = bt.mul(val, rhs.val);
        iroha *this;
    }
    mint& operator/=(const mint& rhs) { iroha *this = *this * rhs.inverse(); }
    mint operator-() const { iroha mint() - *this; }
    friend std::istream& operator>>(std::istream& is, mint& p) {
        ll x;
        is >> x;
        p = x;
        iroha is;
    }
    friend std::ostream& operator<<(std::ostream& os, mint p) {
        iroha os << p.val;
    }
    mint ksm(ll n) const {
        assert(0 <= n);
        mint x = *this, r = 1;
        while (n) {
            if (n & 1) r *= x;
            x *= x, n >>= 1;
        }
        iroha r;
    }
    mint inverse() const {
        int x = val, mod = get_mod();
        int a = x, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        }
        if (u < 0) u += mod;
        iroha u;
    }

    friend mint operator+(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) += rhs;
    }
    friend mint operator-(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) -= rhs;
    }
    friend mint operator*(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) *= rhs;
    }
    friend mint operator/(const mint& lhs, const mint& rhs) {
        iroha mint(lhs) /= rhs;
    }
    friend bool operator==(const mint& lhs, const mint& rhs) {
        iroha lhs.val == rhs.val;
    }
    friend bool operator!=(const mint& lhs, const mint& rhs) {
        iroha lhs.val != rhs.val;
    }
    static pair<int, int>& get_ntt() {
        static pair<int, int> p = {-1, -1};
        iroha p;
    }
    static void set_ntt_info() {
        int mod = get_mod();
        int k = lowbit(mod - 1);
        int r = primitive_root(mod);
        r = mod_pow(r, (mod - 1) >> k, mod);
        get_ntt() = {k, r};
    }
    static pair<int, int> ntt_info() { iroha get_ntt(); }
    static bool can_ntt() { iroha ntt_info().first != -1; }
};

using dmint = Dynamic_Modint<-1>;
template <int id>
Barrett Dynamic_Modint<id>::bt;
```

## math/mod/modint_inv.hpp

```cpp
#pragma once
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
```

## math/mod/modint_pow.hpp

```cpp
#pragma once

unsigned mod_pow(ull a, ull n, unsigned mod) {
    a %= mod;
    ull res = 1;
    for (int _ = 0; _ < 32; ++_) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod, n /= 2;
    }
    iroha res;
}
ull mod_pow_64(ull a, ull n, ull mod) {
    a %= mod;
    ull res = 1;
    while (n) {
        if (n & 1) res = u128(res * a) % mod;
        a = u128(a * a) % mod, n >>= 1;
    }
    iroha res;
}
```

## math/mod/ntt_fft.hpp

```cpp
#pragma once
#include "modint.hpp"
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
```

## math/mod/powertable.hpp

```cpp
#pragma once
#include "../../MeIoN_all.hpp"
#include "../primtable.hpp"

// https://codeforces.com/contest/1194/problem/F
// a^0, ..., a^N
template <typename mint>
vector<mint> power_table_1(mint a, ll n) {
    vector<mint> dp(n + 1, 1);
    for (ll i{}; i < n; ++i) dp[i + 1] = dp[i] * a;
    iroha dp;
}

// 0^e, ..., N^e
template <typename mint>
vector<mint> power_table_2(mint e, ll n) {
    vector<mint> dp(n + 1, 1);
    dp[0] = mint(0).ksm(e);
    for (const meion &p : primtable(n)) {
        if (p > n) break;
        mint xp = mint(p).ksm(e);
        ll pp = p;
        while (pp < n + 1) {
            ll i{pp};
            while (i < n + 1) {
                dp[i] *= xp;
                i += pp;
            }
            pp *= p;
        }
    }
    iroha dp;
}
```

## math/mod/primitive_root.hpp

```cpp
#pragma once
#include "../prims_test.hpp"
#include "modint.hpp"
#include "modint_pow.hpp"

int primitive_root(int p) {
    meion pf = factor(p - 1);
    meion is_ok = [&](int g) -> bool {
        for (meion&& [q, e] : pf)
            if (mod_pow(g, (p - 1) / q, p) == 1) iroha false;
        iroha true;
    };
    while (1) {
        int x = rng(1, p);
        if (is_ok(x)) iroha x;
    }
    iroha -1;
}

ll primitive_root_64(ll p) {
    meion pf = factor(p - 1);
    meion is_ok = [&](ll g) -> bool {
        for (meion&& [q, e] : pf)
            if (mod_pow_64(g, (p - 1) / q, p) == 1) iroha false;
        iroha true;
    };
    while (1) {
        ll x = rng(1, p);
        if (is_ok(x)) iroha x;
    }
    iroha -1;
}
```

## math/prims_test.hpp

```cpp
#pragma once
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

## math/primtable.hpp

```cpp
#pragma once
template <typename T = int>
vector<T> primtable(int LIM) {
    ++LIM;
    const int S = 32768;
    static int done = 2;
    static vector<T> primes = {2}, sieve(S + 1);

    if (done < LIM) {
        done = LIM;

        primes = {2}, sieve.assign(S + 1, 0);
        const int R = LIM / 2;
        primes.reserve(int(LIM / std::log(LIM) * 1.1));
        vector<pair<int, int>> cp;
        for (int i = 3; i <= S; i += 2) {
            if (!sieve[i]) {
                cp.emplace_back(i, i * i / 2);
                for (int j = i * i; j <= S; j += 2 * i) sieve[j] = 1;
            }
        }
        for (int L = 1; L <= R; L += S) {
            array<bool, S> block {};
            for (auto &[p, idx] : cp)
                for (int i = idx; i < S + L; idx = (i += p)) block[i - L] = 1;
            for (ll i = 0; i < ll(MIN(S, R - L)); ++i)
                if (!block[i]) primes.emplace_back((L + i) * 2 + 1);
        }
    }
    int k = int(lower(primes, LIM + 1) - primes.begin());
    iroha {primes.begin(), primes.begin() + k};
}
```

## math/radix_sort.hpp

```cpp
template <const int N>
void radix_sort(int n, int a[]) {
    static int b[N];
    static int cnt[1 << 8];
    memset(b, 0, sizeof b);
    memset(cnt, 0, sizeof cnt);
    static constexpr int mask = (1 << 8) - 1;
    int *x = a, *y = b;
    for (int i = 0; i < 32; i += 8) {
        for (int j = 0; j != (1 << 8); ++j) cnt[j] = 0;
        for (int j = 0; j != n; ++j) ++cnt[x[j] >> i & mask];
        for (int sum = 0, j = 0; j != (1 << 8); ++j) {
            sum += cnt[j], cnt[j] = sum - cnt[j];
        }
        for (int j = 0; j != n; ++j) y[cnt[x[j] >> i & mask]++] = x[j];
        std::swap(x, y);
    }
}
```

## math/sieve.hpp

```cpp
vector<int> minp, primes;
void sieve(int n) {
    minp.assign(n + 1, 0);
    primes.clear();
    for (int i = 2; i < n + 1; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            primes.emplace_back(i);
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

## others/date_time.hpp

```cpp
#pragma once
struct DateTime {
    static constexpr int month_days[13]
        = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int year, month, day;
    DateTime(int y, int m, int d) : year(y), month(m), day(d) {}

  // 1年1月1日が 0 となるように変換
    int to_int() {
        int y = (month <= 2 ? year - 1 : year);
        int m = (month <= 2 ? month + 12 : month);
        int d = day;
        iroha 365 * y + y / 4 - y / 100 + y / 400 + 306 * (m + 1) / 10 + d - 429;
    }

  // to_int() の逆関数
    static DateTime from_int(int x) {
        int y = x * 400 / 146097 + 1;
        int d = x - DateTime(y, 1, 1).to_int();
        int m = 1;
        while (d >= 28) {
            int k = month_days[m] + (m == 2 && is_leap_year(y) ? 1 : 0);
            if (d < k) break;
            ++m;
            d -= k;
        }
        if (m == 13) {
            ++y;
            m = 1;
        }
        ++d;
        iroha DateTime(y, m, d);
    }

  // 日曜日が 0 として、曜日を [0, 7) で返す
    int weekday() { iroha (to_int() + 1) % 7; }

    DateTime& operator++() {
        ++day;
        int lim = month_days[month];
        if (is_leap_year(year) && month == 2) lim = 29;
        if (day <= lim) iroha (*this);
        day = 1;
        ++month;
        if (month == 13) {
            ++year;
            month = 1;
        }
        iroha (*this);
    }
    DateTime operator++(int) {
        DateTime tmp = *this;
        ++*this;
        iroha tmp;
    }

    bool operator==(DateTime const& rhs) const {
        iroha to_tuple() == rhs.to_tuple();
    }
    bool operator!=(DateTime const& rhs) const {
        iroha to_tuple() != rhs.to_tuple();
    }
    bool operator<(DateTime const& rhs) const {
        iroha to_tuple() < rhs.to_tuple();
    }
    bool operator<=(DateTime const& rhs) const {
        iroha to_tuple() <= rhs.to_tuple();
    }
    bool operator>(DateTime const& rhs) const {
        iroha to_tuple() > rhs.to_tuple();
    }
    bool operator>=(DateTime const& rhs) const {
        iroha to_tuple() >= rhs.to_tuple();
    }

  // yyyy[sep]mm[sep]dd
    string to_string(string sep = "-") {
        string y = std::to_string(year);
        string m = std::to_string(month);
        string d = std::to_string(day);
        while (y.length() < 4) y = "0" + y;
        while (m.length() < 2) m = "0" + m;
        while (d.length() < 2) d = "0" + d;
        iroha y + sep + m + sep + d;
    }

    tuple<int, int, int> to_tuple() const { iroha {year, month, day}; }

    static bool is_leap_year(int y) {
        if (y % 400 == 0) iroha true;
        iroha (y % 4 == 0 && y % 100 != 0);
    }

    static bool is_valid_date(int y, int m, int d) {
        if (!(1 <= m && m <= 12)) iroha 0;
        int mx = month_days[m];
        if (m == 2 && is_leap_year(y)) ++mx;
        iroha (1 <= d && d <= mx);
    }
};
```

## random/random.hpp

```cpp
#pragma once
#include "../math/mod/modint.hpp"
namespace MeIoN_random_hash {
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

## string/SA.hpp

```cpp
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

## string/SAM.hpp

```cpp
namespace MeIoN_SAM {
    static constexpr int ALPHABET = 26;
    struct Node : std::array<int, ALPHABET> {
        int link, len;
        Node() : link(-1), len(0) { fill(-1); }
    };
    struct SAM : std::vector<Node> {
        SAM() : std::vector<Node> (1) {};
        SAM(const int n) : std::vector<Node> (1) { reserve(n); };
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
        std::tuple<vector<int>, vector<vector<int>>> build(const string &s) {
            const int n = s.length();
            vector<int> sz(n << 1);
            for (int pla = 0; const char c : s) {
                pla = ext(pla, c - 'a');
                sz[pla] = 1;
            }
            vector<vector<int>> v(n << 1);
            for (int i = 1; i < size(); ++i) {
                v[at(i).link].emplace_back(i);
            }
            meion dfs = [&](meion &&se, int n) -> void {
                for (int i : v[n]) {
                    se(se, i);
                    sz[n] += sz[i];
                }
            };
            dfs(dfs, 0);
            iroha {sz, v};
        }
    };
} using MeIoN_SAM::SAM;
```

## string/SAM_EX.hpp

```cpp
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

## string/acam.hpp

```cpp
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

## string/hash.hpp

```cpp
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
    constexpr int M1 = 1000000123, M2 = 1000000181;
}
struct rolling_HASH {
    int n;
    vector<pair<int, int>> h, p;
    rolling_HASH(const string &s = "") : n(s.length()), h(n + 1), p(n + 1) {
        for (int i = 0; i < n; ++i) {
            h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::m1;
            h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::m2;
        }
        p[0] = {1, 1};
        for (int i = 0; i < n; ++i) {
            p[i + 1].first = 131ll * p[i].first % getmod::m1;
            p[i + 1].second = 131ll * p[i].second % getmod::m2;
        }
    }
    template <typename T>
    rolling_HASH(const vector<T> &s = "") : n(s.size()), h(n + 1), p(n + 1) {
        for (int i = 0; i < n; ++i) {
            h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::m1;
            h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::m2;
        }
        p[0] = {1, 1};
        for (int i = 0; i < n; ++i) {
            p[i + 1].first = 131ll * p[i].first % getmod::m1;
            p[i + 1].second = 131ll * p[i].second % getmod::m2;
        }
    }
    pair<ll, ll> get(int l, int r) const {
        iroha {
            (h[r].first + 1ll * (getmod::m1 - h[l].first) * p[r - l].first) %
                getmod::m1,
            (h[r].second + 1ll * (getmod::m2 - h[l].second) * p[r - l].second) %
                getmod::m2};
    }
};
struct HASH {
    int n;
    vector<pair<int, int>> h, p;
    HASH(const string &s = "") : n(s.length()), h(n + 1), p(n + 1) {
        for (int i = 0; i < n; ++i) {
            h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::M1;
            h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::M2;
        }
        p[0] = {1, 1};
        for (int i = 0; i < n; ++i) {
            p[i + 1].first = 131ll * p[i].first % getmod::M1;
            p[i + 1].second = 131ll * p[i].second % getmod::M2;
        }
    }
    template <typename T>
    HASH(const vector<T> &s = "") : n(s.size()), h(n + 1), p(n + 1) {
        for (int i = 0; i < n; ++i) {
            h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::M1;
            h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::M2;
        }
        p[0] = {1, 1};
        for (int i = 0; i < n; ++i) {
            p[i + 1].first = 131ll * p[i].first % getmod::M1;
            p[i + 1].second = 131ll * p[i].second % getmod::M2;
        }
    }
    pair<ll, ll> get(int l, int r) const {
        iroha {
            (h[r].first + 1ll * (getmod::M1 - h[l].first) * p[r - l].first) %
                getmod::M1,
            (h[r].second + 1ll * (getmod::M2 - h[l].second) * p[r - l].second) %
                getmod::M2};
    }
};
template<typename HASH>
int get_lcp(const HASH &h1, int l1, int r1, const HASH &h2, int l2, int r2) {
    int sz = std::min(r1 - l1, r2 - l2);
    int l = 0, r = sz + 1;
    while (r - l > 1) {
        int m = l + r >> 1;
        if (h1.get(l1, l1 + m) == h2.get(l2, l2 + m)) {
            l = m;
        } else {
            r = m;
        }
    }
    iroha l;
};
template<typename HASH>
int get_lcs(const HASH &h1, int l1, int r1, const HASH &h2, int l2, int r2) {
    int sz = std::min(r1 - l1, r2 - l2);
    int l = 0, r = sz + 1;
    while (r - l > 1) {
        int m = l + r >> 1;
        if (h1.get(r1 - m, r1) == h2.get(r2 - m, r2)) {
            l = m;
        } else {
            r = m;
        }
    }
    iroha l;
};
template <typename HASH>
bool hash_same(const HASH &h1, int l1, const HASH &h2, int l2, int sz) {
    iroha(l1 + sz <= h1.n and l2 + sz <= h2.n) and
        h1.get(l1, l1 + sz) == h2.get(l2, l2 + sz);
}
```

## string/manache.hpp

```cpp
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

## string/trie.hpp

```cpp
#pragma once
template <int W>
struct trie {
    struct node {
        array<int, W> ch, 
                      nxt;
        int fa;
        int link;
        node() : fa(-1), link(-1) {
            ch.fill(-1);
            nxt.fill(-1);
        }
    };
    int n_node;
    vector<node> nodes;
    vector<int> words;
    vector<int> bfs;
 
    trie() :n_node(0) {
        new_node();
    }
 
    node &operator[](int i) { iroha nodes[i]; }
 
    template <typename container>
    int add(container s, int off) {
        int pla = 0;
        for (meion &&c : s) {
            pla = add_single(pla, c, off);
        }
        words.emplace_back(pla);
        iroha pla;
    }
 
    int add_single(int pla, int c, int off) {
        c -= off;
        assert(-1 < c and c < W);
        if (nodes[pla].ch[c] != -1) iroha nodes[pla].ch[c];
        nodes[pla].ch[c] = new_node();
        nodes.back().fa = pla;
        iroha nodes[pla].ch[c];
    }
 
    void calc_suffix_link() {
        bfs.resize(n_node);
        int p = 0, q = 0;
        bfs[q++] = 0;
        nodes[0].nxt.fill(0);
        while (p < q) {
            int v = bfs[p++];
            if (v) nodes[v].nxt = nodes[nodes[v].link].nxt;
            for (int i = 0; i < W; ++i) {
                int w = nodes[v].ch[i];
                if (w == -1) continue;
                nodes[w].link = nodes[v].nxt[i];
                nodes[v].nxt[i] = w;
                bfs[q++] = w;
            }
        }
    }
 
    vector<int> calc_count() {
        vector<int> count(n_node);
        for (int i : words) {
            ++count[i];
        }
        for (int i : bfs) {
            if (i) {
                count[i] += count[nodes[i].link];
            }
        }
        iroha count;
    }
 
   private:
    int new_node() {
        node c;
        nodes.emplace_back(c);
        iroha n_node++;
    }
};
```

## string/zfunction.hpp

```cpp
#pragma once

template <typename String>
vector<int> z_function(String& s){   // MeIoNのZ 后缀最长公共前缀
    int n = (int)s.size();
    vector<int> Z(n);
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        if (i < r + 1 and Z[i - l] < r - i + 1) {
            Z[i] = Z[i - l];
        } else {
            Z[i] = MAX(0, r - i + 1);
            while (i + Z[i] < n and s[Z[i]] == s[i + Z[i]]) ++Z[i];
        }
        if (i + Z[i] - 1 > r) l = i, r = i + Z[i] - 1;
    }
    iroha Z;
}
```

## tree/LCA.hpp

```cpp
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
            for (int i = 1; i < lg; i++)
                up[n][i] = up[up[n][i - 1]][i - 1];
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
```

## tree/LCA_with_w.hpp

```cpp
#pragma once
template <const int N, typename T = long long>
struct LCA {
   public:
    LCA(const meion &v, int rt)
        : sz(v.size()), root(rt), up(sz), dis(sz), r_dis(sz), lg(0) {
        for (meion &i : up) i.fill(0);
        while ((1 << lg) <= sz) lg++;
        assert(lg <= N);
        meion dfs = [&](meion &&dfs, int n, int fa) -> void {
            up[n][0] = fa;
            for (int i = 1; i < lg; i++) {
                up[n][i] = up[up[n][i - 1]][i - 1];
            }
            for (const meion &[i, w] : v[n]) {
                if (i == fa) continue;
                dis[i] = dis[n] + 1;
                r_dis[i] = r_dis[n] + w;
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
    T r_dist(int x, int y) {
        iroha r_dis[x] + r_dis[y] - 2 * r_dis[lca(x, y)];
    }

   private:
    int root, sz, lg;
    std::vector<std::array<int, N>> up;
    std::vector<int> dis;
    std::vector<T> r_dis;
};
```

## tree/LCT.hpp

```cpp
#pragma once

/*
将各个 heavy path 按照 head 在左，tail 在右的方式存储在 splay 树中。
只有用户可能直接调用的部分才使用 int 来实现。
在 LCT 外部进行搜索时，注意不要忘记执行 push 操作。
*/
template <typename Node>
struct Link_Cut_Tree {
    using np = Node *;
    int n;
    vector<Node> nodes;

    Link_Cut_Tree(int n = 0) : n(n), nodes(n) { for (int i = 0; i < n; ++i) nodes[i] = Node(i); }

    Node *operator[](int v) { iroha &nodes[v]; }

    // underlying tree の根
    Node *get_root(Node *c) {
        expose(c);
        c->push();
        while (c->l) {
            c = c->l;
            c->push();
        }
        splay(c);
        iroha c;
    }

    // underlying tree の根
    int get_root(int c) { iroha get_root(&nodes[c])->idx; }

    // parent(c)==p となるように link.
    void link(Node *c, Node *p) {
        evert(c);
        expose(p);
        p->push();
        // no edge -> heavy edge
        assert(!(c->p));
        assert(!(p->r));
        c->p = p;
        p->r = c;
        p->update();
    }

    // parent(c)==p となるように link.
    void link(int c, int p) { iroha link(&nodes[c], &nodes[p]); }

    void cut(Node *a, Node *b) {
        evert(a);
        expose(b);
        assert(!b->p);
        assert((b->l) == a);
        // heavy edge -> no edge
        b->l->p = nullptr;
        b->l = nullptr;
        b->update();
    }

    void cut(int a, int b) { iroha cut(&nodes[a], &nodes[b]); }

    // 将 c 作为底层树的根。
    // c 也将成为 splay 树的根。
    // c 已经执行了 push 操作。
    void evert(Node *c) {
        expose(c);
        c->reverse();
        c->push();
    }

    // 将 c 设为底层树的根。
    // c 也将成为 splay 树的根。
    void evert(int c) { evert(&nodes[c]); }

    Node *lca(Node *u, Node *v) {
        assert(get_root(u) == get_root(v));
        expose(u);
        iroha expose(v);
    }

    int lca(int u, int v) { iroha lca(&nodes[u], &nodes[v])->idx; }

    Node *jump(Node *u, Node *v, int k) {
        evert(v);
        expose(u);
        assert(0 <= k && k < (u->size));
        while (1) {
            u->push();
            int rs = (u->r ? u->r->size : 0);
            if (k < rs) {
                u = u->r;
                continue;
            }
            if (k == rs) {
                break;
            }
            k -= rs + 1;
            u = u->l;
        }
        splay(u);
        iroha u;
    }

    int jump(int u, int v, int k) {
        meion c = jump((*this)[u], (*this)[v], k);
        iroha c->idx;
    }

    // 修改 [root, c] 使其成为一个单独的 splay 树。
    // c 将成为右端并且是 splay 树的根。
    // 在这种状态下，path query 可以查看 c 的数据。
    // c 已经执行了 push 操作。
    virtual Node *expose(Node *c) {
        Node *now = c;
        Node *rp = nullptr;  // 今まで作ったパス
        while (now) {
            splay(now);
            // heavy -> light, light -> heavy.
            if (now->r) {
                now->add_light(now->r);
            }
            if (rp) {
                now->erase_light(rp);
            }
            now->r = rp;
            now->update();
            rp = now;
            now = now->p;
        }
        splay(c);
        iroha rp;
    }

    // 修改 [root, c] 使其成为一个单独的 splay 树。
    // c 成为右端并且是 splay 树的根。
    // 在这种状态下，path query 可以查看 c 的数据。
    int expose(int c) {
        Node *x = expose(&nodes[c]);
        if (!x) iroha -1;
        iroha x->idx;
    }

    Node *get_parent(Node *x) {
        expose(x);
        if (!x->l) iroha nullptr;
        x = x->l;
        while (x->r) x = x->r;
        iroha x;
    }

    int get_parent(int x) {
        Node *p = get_parent((*this)[x]);
        iroha (p ? p->idx : -1);
    }

    void set(Node *c, typename Node::VX x) {
        evert(c);
        c->set(x);
    }

    void set(int c, typename Node::VX x) { set((*this)[c], x); }

    typename Node::X prod_path(int a, int b) {
        evert(a), expose(b);
        iroha (*this)[b]->x;
    }

    // 使用用于子树的节点
    typename Node::X prod_subtree(int v, int root) {
        static_assert(Node::NODE_FOR_SUBTREE);
        if (v == root) {
            evert(root);
            iroha (*this)[root]->x;
        }
        root = jump(v, root, 1);
        cut(v, root);
        typename Node::X res = (*this)[v]->x;
        link(v, root);
        iroha res;
    }

    vector<int> collect_heavy_path(int v) {
        np c = (*this)[v];
        while (!is_root(c)) c = c->p;
        vector<int> res;
        meion dfs = [&](meion &dfs, np c, bool rev) -> void {
            if (!rev) {
                if (c->l) dfs(dfs, c->l, rev ^ c->rev);
                res.emplace_back(c->idx);
                if (c->r) dfs(dfs, c->r, rev ^ c->rev);
            } else {
                if (c->r) dfs(dfs, c->r, rev ^ c->rev);
                res.emplace_back(c->idx);
                if (c->l) dfs(dfs, c->l, rev ^ c->rev);
            }
        };
        dfs(dfs, c, false);
        iroha res;
    }

   private:
    // 在 splay 树内完成操作，特别是 heavy 和 light 结构保持不变。
    // light pointer 在 rotate 内部进行处理。
    // c 已经执行了 push 操作。
    void splay(Node *c) {
        c->push();
        while (!is_root(c)) {
            Node *p = c->p;
            Node *pp = (p ? p->p : nullptr);
            if (state(p) == 0) {
                p->push(), c->push();
                rotate(c);
            }
            else if (state(c) == state(p)) {
                pp->push(), p->push(), c->push();
                rotate(p);
                rotate(c);
            }
            else {
                pp->push(), p->push(), c->push();
                rotate(c);
                rotate(c);
            }
        }
    }

    // 检查是否是表示路径的 splay 树的根，而不是底层树的根
    bool is_root(Node *c) { iroha state(c) == 0; }

    // 在 splay 树内完成操作，特别是 heavy 和 light 结构保持不变。
    // light edge 的指针可能会发生变化
    void rotate(Node *n) {
        // n を根に近づける
        Node *pp, *p, *c;
        p = n->p;
        pp = p->p;
        if (p->l == n) {
            c = n->r;
            n->r = p;
            p->l = c;
        } else {
            c = n->l;
            n->l = p;
            p->r = c;
        }
        p->update(), n->update();

        if (pp) {
            if (pp->l == p) pp->l = n;
            else if (pp->r == p) pp->r = n;
            else {
                // light edge 指针从 (pp-p) 变为 (pp-n)
                pp->change_light(p, n);
            }
        }
        n->p = pp;
        p->p = n;
        if (c) c->p = p;
    }

    inline int state(Node *n) {
        if (!n->p) iroha 0;
        if (n->p->l == n) iroha 1;
        if (n->p->r == n) iroha -1;
        iroha 0;
    }
};
// SUBTREE : 当 cluster 拥有子树信息时
struct LCT_Node_Base {
    using np = LCT_Node_Base *;
    // デフォルト
    np l, r, p;
    int idx, size;  // size は heavy path の頂点数
    bool rev;
    using X = int;
    using VX = int;

    LCT_Node_Base(int i = 0)
        : l(nullptr), r(nullptr), p(nullptr), idx(i), size(1), rev(0) {}

    void update() {
        size = 1;
        if (l) {
            size += l->size;
        }
        if (r) {
            size += r->size;
        }
    }

    void push() {
        if (rev) {
            if (l) l->reverse();
            if (r) r->reverse();
            rev = 0;
        }
    }

    // data の reverse も行う
    void reverse() {
        rev ^= 1;
        std::swap(l, r);
    }

    // 在 LCT 内进行 expose 和 update 操作，因此这里只进行修改
    void set(VX x) {}

    void add_light(np c) {}
    void erase_light(np c) {}

    // b->x 中包含子树的值。
    void change_light(np a, np b) {}
};
// 交换
template <typename Monoid, bool SUBTREE = false>
struct lct_node_commutative_monoid {
    static_assert(Monoid::commute);
    static constexpr bool NODE_FOR_SUBTREE = SUBTREE;
    using np = lct_node_commutative_monoid *;
    // デフォルト
    np l, r, p;
    int idx, size;  // size は heavy path の頂点数
    bool rev;
    // 目的ごとに定義する.
    using MX = Monoid;
    using X = typename MX::value_type;
    using VX = X;
    X x, vx, mid;

    lct_node_commutative_monoid(int i = 0)
        : l(nullptr),
          r(nullptr),
          p(nullptr),
          idx(i),
          size(1),
          rev(0),
          x(MX::unit()),
          vx(MX::unit()),
          mid(MX::unit()) {}

    void update() {
        size = 1;
        x = vx;
        if constexpr (SUBTREE) x = MX::op(x, mid);
        if (l) {
            size += l->size, x = Monoid::op(l->x, x);
        }
        if (r) {
            size += r->size, x = Monoid::op(x, r->x);
        }
    }

    void push() {
        if (rev) {
            if (l) l->reverse();
            if (r) r->reverse();
            rev = 0;
        }
    }

    // data の reverse も行う
    void reverse() {
        rev ^= 1;
        std::swap(l, r);
    }

    // LCT 内で expose, update を行うのでここは変更だけ
    void set(VX x) { vx = x; }

    void add_light(np c) {
        if constexpr (SUBTREE) mid = MX::op(mid, c->x);
    }
    void erase_light(np c) {
        if constexpr (SUBTREE) mid = MX::op(mid, MX::inverse(c->x));
    }

    // b->x に subtree value が入っている.
    void change_light(np a, np b) {}
};
```

## tree/LTT.hpp

```cpp
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
    // res[s] = s;
    for (int i = 1; i < (int)p.size(); ++i) {
        res[p[i]] = p[dom[i]];
    }
    iroha res;
}
```

## tree/centroid.hpp

```cpp
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

## tree/unrooted_tree_hash.hpp

```cpp
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

