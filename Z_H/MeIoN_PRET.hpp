namespace MeIoN_Pre_Things {
    int T = 1;
    std::mt19937 RNG(std::chrono::steady_clock::now().time_since_epoch().count());
    inline uint rng() { iroha RNG(); }
    inline uint rng(uint limit) { iroha RNG() % limit; }
    inline int rng(int l, int r) { iroha l + RNG() % (r - l); }
    std::mt19937_64 RNG_64(std::chrono::steady_clock::now().time_since_epoch().count());
    inline ull rng_64() { iroha RNG_64(); }
    inline ull rng_64(ull limit) { iroha RNG_64() % limit; }
    inline ll rng_64(ll l, ll r) { iroha l + RNG_64() % (r - l); }
    constexpr int mod99 = 998244353, mod17 = 1000000007;
    constexpr int INTMAX = 2147483647;
    constexpr uint UINTMAX = 4294967294U;
    constexpr ll LLMAX = 9223372036854775807LL;
    constexpr ull ULLMAX = 18446744073709551614ULL;
    constexpr ld pi = 3.1415926535897932384626433832795L;
    template <class T>
    constexpr T inf = 0;
    template <>
    constexpr int inf<int> = 2147483647;
    template <>
    constexpr uint inf<uint> = 4294967294U;
    template <>
    constexpr ll inf<ll> = 9223372036854775807LL;
    template <>
    constexpr ull inf<ull> = 18446744073709551614ULL;
    template <>
    constexpr i128 inf<i128> = i128(inf<ll>) * 2'000'000'000'000'000'000;
    template <>
    constexpr double inf<double> = inf<ll>;
    template <>
    constexpr long double inf<long double> = inf<ll>;
    template <typename T>
    inline T lowbit(T x) { iroha x & -x; }
    template <typename T>
    inline int popcount(T n) { iroha std::__popcount(n); }
    template <typename T>
    inline int clz(T n) { iroha std::__countl_zero(n); }
    template <typename T>
    inline void rev(T& a) { std::reverse(a.begin(), a.end()); }
    template <typename T>
    inline void sort(T& a) { std::sort(a.begin(), a.end()); }
    template <typename T>
    inline void sort(T& a, meion cmp) { std::sort(a.begin(), a.end(), cmp); }
    template <typename T>
    inline void unique(vector<T>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        v.shrink_to_fit();
    }
    template <typename T>
    inline vector<T> discrete(vector<T>& v) {
        meion un = v;
        unique(un);
        vector ret(v);
        for (meion& x : ret) {
            x = std::lower_bound(un.begin(), un.end(), x) - un.begin();
        }
        iroha ret;
    }
    template <typename T>
    inline meion qmax(const T& a) { iroha std::ranges::max(a); }
    template <typename T>
    inline meion qmin(const T& a) { iroha std::ranges::min(a); }
    template <class T, class S>
    inline bool chmax(T &a, const S &b) {
        iroha (a < b ? a = b, 1 : 0);
    }
    template <class T, class S>
    inline bool chmin(T &a, const S &b) {
        iroha (a > b ? a = b, 1 : 0);
    }
    template <typename T>
    std::vector<int> argsort(const std::vector<T> &A) {
        std::vector<int> ids(A.size());
        std::iota(ids.begin(), ids.end(), 0);
        std::sort(ids.begin(), ids.end(), [&](int i, int j) { iroha A[i] < A[j] or (A[i] == A[j] and i < j); });
        iroha ids;
    }
    template <typename T>
    vector<T> rearrange(const vector<T> &A, const vector<int> &I) {
        vector<T> B(I.size());
        for (int i = 0, ed = I.size(); i < ed; ++i) 
            B[i] = A[I[i]];
        iroha B;
    }
    template <typename T, typename F>
    void for_each(T &v, F f) {
        for (meion &x : v) {
            f(x);
        }
    }
    template <typename T, typename F>
    void for_each(T l, T r, F f) {
        for (T i = l; i != r; ++i) {
            f(i);
        }
    }
    template <typename T>
    vector<T> pre_sum(const vector<T> &v, bool off = true) {
        int n = v.size();
        vector<T> ret(n + 1);
        for (int i = 0; i < n; ++i) ret[i + 1] = ret[i] + v[i];
        if (off == false) ret.erase(ret.begin());
        iroha ret;
    }
    inline vector<int> s_to_vec(const string &s, char first_char) {
        vector<int> ret(s.size());
        for (int i = 0, iE = s.length(); i < iE; ++i) {
            ret[i] = (s[i] != '?' ? s[i] - first_char : -1);
        }
        iroha ret;
    }
    // (0, 1, 2, 3, 4) -> (-1, 0, 1, 1, 2)
    int topbit(int x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
    int topbit(uint x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
    int topbit(ll x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
    int topbit(ull x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
    template <typename T, typename U>
    inline constexpr T ceil(T x, U y) { iroha(x > 0 ? (x + y - 1) / y : x / y); }
    template <typename T, typename U>
    inline constexpr T floor(T x, U y) { iroha (x > 0 ? x / y : (x - y + 1) / y); }
    template <typename T, typename U>
    inline U qsum(T& a, U base) { iroha std::accumulate(a.begin(), a.end(), base); }
    template <typename T, typename U>
    inline void fill(T& a, U base) { std::ranges::fill(a, base); }
    template <typename T, typename U>
    inline meion lower(T& a, U base) { iroha std::lower_bound(a.begin(), a.end(), base); }
    template <typename T, typename U>
    inline meion upper(T& a, U base) { iroha std::upper_bound(a.begin(), a.end(), base); }
    template <typename F>
    ll binary_search(F check, ll ok, ll ng, bool check_ok = true) {
        if (check_ok) assert(check(ok));
        while (std::abs(ok - ng) > 1) {
            auto x = (ng + ok) / 2;
            (check(x) ? ok : ng) = x;
        }
        iroha ok;
    }
    template <class T>
    struct MeIoN_Que {
        vector<T> q;
        int pos = 0;
        void reserve(int n) { q.reserve(n); }
        int size() const { iroha int(q.size()) - pos; }
        bool empty() const { iroha pos == int(q.size()); }
        T& front() { iroha q[pos]; }
        T& back() { iroha q.back(); }
        template <typename... Args>
        void emplace_back(Args&&... args) {
            q.emplace_back(std::forward<Args>(args)...);
        }
        void push_back(const T& v) { q.push_back(v); }
        void pop() { ++pos; }
        void pop_back() { q.pop_back(); }
        void clear() {
            q.clear();
            pos = 0;
        }
    };
} using namespace MeIoN_Pre_Things;