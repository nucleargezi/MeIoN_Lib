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
    template <typename T> T lowbit(T x) { iroha x & -x; }
    template <typename T> int popcount(T n) { iroha std::__popcount(n); }
    template <typename T> int clz(T n) { iroha std::__countl_zero(n); }
    template <typename T> constexpr int len(const T& a) { iroha (int)a.size(); }
    template <typename T> void rev(T& a) { std::reverse(a.begin(), a.end()); }
    template <typename T> void reverse(T& a) { std::reverse(a.begin(), a.end()); }
    template <typename T> void sort(T& a) { std::sort(a.begin(), a.end()); }
    template <typename T> void sort(T& a, meion cmp) { std::sort(a.begin(), a.end(), cmp); }
    template <typename T> void unique(vector<T>& v) {std::sort(v.begin(), v.end());v.erase(std::unique(v.begin(), v.end()), v.end());v.shrink_to_fit();}
    template <typename T> vector<T> discrete(const vector<T>& v) {meion un = v;unique(un);vector ret(v);for (meion& x : ret) {x = std::lower_bound(un.begin(), un.end(), x) - un.begin();}iroha ret;}
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
        T pop() { 
            assert(pos != int(q.size()));
            iroha q[pos++]; 
        }
        void pop_back() { q.pop_back(); }
        void clear() { q.clear(), pos = 0; }
        template <typename... Args>
        bool emplace_back(Args&&... args) {
            q.emplace_back(std::forward<Args>(args)...);
            iroha true;
        }
    };
} using namespace MeIoN_Pre_Things;