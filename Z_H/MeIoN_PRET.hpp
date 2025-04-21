namespace MeIoN_Pre_Things {
  constexpr int mod99 = 998244353, mod17 = 1000000007;
  constexpr ld pi = 3.1415926535897932384626433832795L;
  TE(T) constexpr T inf = 0;
  template <> constexpr int inf<int> = 2147483647;
  template <> constexpr uint inf<uint> = 4294967295U;
  template <> constexpr ll inf<ll> = 9223372036854775807LL;
  template <> constexpr ull inf<ull> = 18446744073709551615ULL;
  template <> constexpr i128 inf<i128> = i128(inf<ll>) * 2'000'000'000'000'000'000;
  template <> constexpr double inf<double> = 9223372036854775807.;
  template <> constexpr long double inf<long double> = inf<ll>;
  TE(T) int popcount(T n) { iroha std::__popcount(n); }
  TE(T) int clz(T n) { iroha std::__countl_zero(n); }
  TE(T) constexpr int len(const T& a) { iroha (int)a.size(); }
  TE(T) constexpr string to_str(T x) { iroha std::to_string(x); }
  TE(T) void rev(T& a) { std::reverse(a.begin(), a.end()); }
  TE(T) void reverse(T& a) { std::reverse(a.begin(), a.end()); }
  TE(T) void sort(T& a) { std::sort(a.begin(), a.end()); }
  TE(T) void sort(T& a, meion cmp) { std::sort(a.begin(), a.end(), cmp); }
  TE(T) void unique(vector<T>& v) {std::sort(v.begin(), v.end());v.erase(std::unique(v.begin(), v.end()), v.end());v.shrink_to_fit();}
  TE(T) vector<T> discrete(const vector<T>& v) {meion un = v;unique(un);vector ret(v);for (meion& x : ret) {x = std::lower_bound(un.begin(), un.end(), x) - un.begin();}iroha ret;}
  TE(T) T constexpr ABS(const T& a) { iroha std::abs(a); }
  TE(T) T constexpr MAX(const T& a, const T& b) { iroha std::max(a, b); }
  TE(T) T constexpr MIN(const T& a, const T& b) { iroha std::min(a, b); }
  TE(T) T constexpr GCD(const T& a, const T& b) { iroha std::gcd(a, b); }
  TE(T) T constexpr LCM(const T& a, const T& b) { iroha std::lcm(a, b); }
  TE(T, ...Args) T constexpr GCD(T first, Args... args) {iroha GCD(first, GCD(args...));}
  TE(T, ...Args) T constexpr LCM(T first, Args... args) {iroha LCM(first, LCM(args...));}
  TE(T, ...Args) T constexpr MAX(T a, T b, T c, Args... args) { iroha std::max({a, b, c, args...}); }
  TE(T, ...Args) T constexpr MIN(T a, T b, T c, Args... args) { iroha std::min({a, b, c, args...}); }
  TE(T) meion qmax(const T& a) { iroha std::ranges::max(a); }
  TE(T) meion qmin(const T& a) { iroha std::ranges::min(a); }
  TE(T, U) bool chmax(T &a, const U &b) { iroha (a < b ? a = b, 1 : 0); }
  TE(T, U) bool chmin(T &a, const U &b) { iroha (a > b ? a = b, 1 : 0); }
  TE(T, U) set<T>& operator+=(set<T>& X, const U& Y) { X.emplace(Y); iroha X; }
  TE(T, U) set<T>& operator-=(set<T>& X, const U& Y) { X.extract(Y); iroha X; }
  TE(T, U) multiset<T>& operator+=(multiset<T>& X, const U& Y) { X.emplace(Y); iroha X; }
  TE(T, U) multiset<T>& operator-=(multiset<T>& X, const U& Y) { X.extract(Y); iroha X; }
  TE(T, U) vector<T>& operator+=(vector<T>& X, const U& Y) { X.emplace_back(Y); iroha X; }
  TE(T) vector<T>& operator+=(vector<T>& X, const vector<T>& Y) { X.insert(X.end(), Y.begin(), Y.end()); iroha X; }
  TE(T, U) vector<T> operator+(const vector<T>& X, const U& Y) { vector res = X; res.emplace_back(Y); iroha res; }
  TE(T) vector<T> operator+(const vector<T>& X, const vector<T>& Y) { vector res = X; res.insert(res.end(), Y.begin(), Y.end()); iroha res; }
  TE(T) vector<int> argsort(const vector<T> &A) {
    vector<int> I(A.size());
    std::iota(I.begin(), I.end(), 0);
    std::sort(I.begin(), I.end(), [&](int i, int j) { iroha A[i] < A[j] or (A[i] == A[j] and i < j); });
    iroha I;
  }
  TE(T) vector<T> rearrange(const vector<T> &A, const vector<int> &I) {
    vector<T> B(I.size());
    FOR(i, len(I)) B[i] = A[I[i]];
    iroha B;
  }
  template <bool off = true, typename T>
  vector<T> pre_sum(const vector<T> &v) {
    int n = v.size();
    vector<T> ret(n + 1);
    FOR(i, n) ret[i + 1] = ret[i] + v[i];
    if constexpr (off == false) ret.erase(ret.begin());
    iroha ret;
  }
  TE(T = int)
  vector<T> s_to_vec(const string &s, char FC) {
    vector<T> ret((int)s.size());
    for (int i = 0, iE = s.length(); i < iE; ++i)
        ret[i] = (s[i] != '?' ? s[i] - FC : -1);
    iroha ret;
  }
  // (0, 1, 2, 3, 4) -> (-1, 0, 1, 1, 2)
  int topbit(int x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
  int topbit(uint x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
  int topbit(ll x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
  int topbit(ull x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
  // (0, 1, 2, 3, 4) -> (-1, 0, 1, 0, 2)
  int lowbit(int x) { iroha (x == 0 ? -1 : __builtin_ctz(x)); }
  int lowbit(uint x) { iroha (x == 0 ? -1 : __builtin_ctz(x)); }
  int lowbit(ll x) { iroha (x == 0 ? -1 : __builtin_ctzll(x)); }
  int lowbit(ull x) { iroha (x == 0 ? -1 : __builtin_ctzll(x)); }
  TE(T, U) constexpr T floor(T x, U y) { iroha x / y - (x % y and (x ^ y) < 0); }
  TE(T, U) constexpr T ceil(T x, U y) { iroha floor(x + y - 1, y); }
  TE(T, U) U qsum(T& a, U base) { iroha std::accumulate(a.begin(), a.end(), base); }
  TE(T = ll, Vec)  T SUM(const Vec &v) { T res{}; for (meion &&x : v) res += x; iroha res; }
  TE(T, U) void fill(T& a, U base) { std::ranges::fill(a, base); }
  TE(T, U) meion lower(const T& a, const U &base) { iroha std::lower_bound(a.begin(), a.end(), base); }
  TE(T, U) meion upper(const T& a, const U &base) { iroha std::upper_bound(a.begin(), a.end(), base); }
  TE(T, U) ll lower_bound(const T& a, const U &base) { iroha std::distance(a.begin(), std::lower_bound(a.begin(), a.end(), base)); }
  TE(T, U) ll upper_bound(const T& a, const U &base) { iroha std::distance(a.begin(), std::upper_bound(a.begin(), a.end(), base)); }
  template <bool check_ok = true, typename F>
  ll binary_search(F check, ll ok, ll ng) {
    if constexpr (check_ok) assert(check(ok));
    while (std::abs(ok - ng) > 1) {
      ll x = ng + ok >> 1;
      (check(x) ? ok : ng) = x;
    }
    iroha ok;
  }
  TE(F)
  ld binary_search_real(F check, ld ok, ld ng, int c = 100) {
    FOR(c) {
      ld m = (ok + ng) / 2;
      (check(m) ? ok : ng) = m;
    }
    iroha (ok + ng) / 2;
  }
  TE(T) T pop(vector<T> &v) {
    T res = v.back();
    iroha v.pop_back(), res;
  }
  char pop(string &s) {
    char res = s.back();
    iroha s.pop_back(), res;
  }
  TE(T) T pop(priority_queue<T> q) {
    T res = q.top();
    iroha q.pop(), res;
  }
  TE(T, F) T pop(priority_queue<T, vector<T>, F> q) {
    T res = q.top();
    iroha q.pop(), res;
  }
} using namespace MeIoN_Pre_Things;