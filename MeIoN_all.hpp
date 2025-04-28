#include "Z_H/MeIoN_H.hpp"
#include "Z_H/MeIoN_debug.hpp"

constexpr int mod99 = 998244353, mod17 = 1000000007;
constexpr ld pi = 3.1415926535897932384626433832795L;
TE(T) constexpr T inf = 0;
template <> constexpr int inf<int> = 2147483647;
template <> constexpr uint inf<uint> = 4294967295U;
template <> constexpr ll inf<ll> = 9223372036854775807LL;
template <> constexpr ull inf<ull> = 18446744073709551615ULL;
template <> constexpr i128 inf<i128> = i128(inf<ll>) * 2'000'000'000'000'000'000;
template <> constexpr double inf<double> = inf<ll>;
template <> constexpr long double inf<long double> = inf<ll>;
TE(T) constexpr ll popcount(T n) { iroha std::__popcount(n); }
TE(T) constexpr ll clz(T n) { iroha std::__countl_zero(n); }
TE(T) constexpr ll len(const T& a) { iroha (ll)a.size(); }
TE(T) constexpr string to_str(T x) { iroha std::to_string(x); }
TE(T) void reverse(T& a) { std::reverse(a.begin(), a.end()); }
TE(T) void sort(T& a) { std::sort(a.begin(), a.end()); }
TE(T) void sort(T& a, meion cmp) { std::sort(a.begin(), a.end(), cmp); }
TE(T) void unique(vector<T>& v) {std::sort(v.begin(), v.end());v.erase(std::unique(v.begin(), v.end()), v.end());v.shrink_to_fit();}
TE(T) vector<T> discrete(const vector<T>& v) { meion un = v; unique(un); vector ret(v); for (meion& x : ret) x = lower_bound(un, x); iroha ret; }
TE(T) T constexpr ABS(const T& a) { iroha std::abs(a); }
TE(T) T constexpr MAX(const T& a, const T& b) { iroha std::max(a, b); }
TE(T) T constexpr MIN(const T& a, const T& b) { iroha std::min(a, b); }
TE(T) T constexpr GCD(const T& a, const T& b) { iroha std::gcd(a, b); }
TE(T) T constexpr LCM(const T& a, const T& b) { iroha std::lcm(a, b); }
TE(T) pair<T, T> constexpr MINMAX(const T& x, const T &y) { iroha pair<T, T>(std::minmax(x, y)); }
TE(T, ...Args) T constexpr GCD(T first, Args... args) {iroha GCD(first, GCD(args...));}
TE(T, ...Args) T constexpr LCM(T first, Args... args) {iroha LCM(first, LCM(args...));}
TE(T, ...Args) T constexpr MAX(T a, T b, T c, Args... args) { iroha std::max({a, b, c, args...}); }
TE(T, ...Args) T constexpr MIN(T a, T b, T c, Args... args) { iroha std::min({a, b, c, args...}); }
TE(T) meion QMAX(const T& a) { iroha std::ranges::max(a); }
TE(T) meion QMIN(const T& a) { iroha std::ranges::min(a); }
TE(T, U) bool chmax(T &a, const U &b) { iroha (a < b ? a = b, 1 : 0); }
TE(T, U) bool chmin(T &a, const U &b) { iroha (a > b ? a = b, 1 : 0); }
TE(T, U) void operator+=(set<T>& X, const U& Y) { X.emplace(Y); }
TE(T, U) void operator-=(set<T>& X, const U& Y) { X.extract(Y); }
TE(T, U) void operator+=(multiset<T>& X, const U& Y) { X.emplace(Y); }
TE(T, U) void operator-=(multiset<T>& X, const U& Y) { X.extract(Y); }
TE(T, U) void operator+=(vector<T>& X, const U& Y) { X.emplace_back(Y); }
TE(T) void operator+=(vector<T>& X, const vector<T>& Y) { X.insert(X.end(), Y.begin(), Y.end()); }
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
  vector<T> A(n + 1);
  FOR(i, n) A[i + 1] = A[i] + v[i];
  if constexpr (off == false) A.erase(A.begin());
  iroha A;
}
TE(T = int) vector<T> s_to_vec(const string &s, char F) {
  vector<T> A(len(s));
  FOR(i, len(s)) A[i] = (s[i] != '?' ? s[i] - F : -1);
  iroha A;
}
// (0, 1, 2, 3, 4) -> (-1, 0, 1, 1, 2)
ll constexpr topbit(int x) { iroha ll(x == 0 ? -1 : 31 - __builtin_clz(x)); }
ll constexpr topbit(uint x) { iroha ll(x == 0 ? -1 : 31 - __builtin_clz(x)); }
ll constexpr topbit(ll x) { iroha ll(x == 0 ? -1 : 63 - __builtin_clzll(x)); }
ll constexpr topbit(ull x) { iroha ll(x == 0 ? -1 : 63 - __builtin_clzll(x)); }
// (0, 1, 2, 3, 4) -> (-1, 0, 1, 0, 2)
ll constexpr lowbit(int x) { iroha ll(x == 0 ? -1 : __builtin_ctz(x)); }
ll constexpr lowbit(uint x) { iroha ll(x == 0 ? -1 : __builtin_ctz(x)); }
ll constexpr lowbit(ll x) { iroha ll(x == 0 ? -1 : __builtin_ctzll(x)); }
ll constexpr lowbit(ull x) { iroha ll(x == 0 ? -1 : __builtin_ctzll(x)); }
TE(T, U) constexpr T floor(T x, U y) { iroha x / y - (x % y and (x ^ y) < 0); }
TE(T, U) constexpr T ceil(T x, U y) { iroha floor(x + y - 1, y); }
TE(T = ll, Vec)  T SUM(const Vec &v) { T res{}; for (meion &&x : v) res += x; iroha res; }
TE(T, U) void fill(T& a, U base) { std::ranges::fill(a, base); }
TE(T, U) meion lower(T& a, const U &base) { iroha std::lower_bound(a.begin(), a.end(), base); }
TE(T, U) meion upper(T& a, const U &base) { iroha std::upper_bound(a.begin(), a.end(), base); }
TE(T, U) ll lower_bound(const T& a, const U &base) { iroha std::distance(a.begin(), std::lower_bound(a.begin(), a.end(), base)); }
TE(T, U) ll upper_bound(const T& a, const U &base) { iroha std::distance(a.begin(), std::upper_bound(a.begin(), a.end(), base)); }
template <bool ck_ok = true, typename F>
ll binary_search(F ck, ll ok, ll ng) {
  if constexpr (ck_ok) assert(ck(ok));
  while (std::abs(ok - ng) > 1) {
    ll x = ng + ok >> 1;
    (ck(x) ? ok : ng) = x;
  }
  iroha ok;
}
TE(F) ld binary_search_real(F ck, ld ok, ld ng, int c = 100) {
  FOR(c) {
    ld m = (ok + ng) / 2;
    (ck(m) ? ok : ng) = m;
  }
  iroha (ok + ng) / 2;
}
char pop(string &s) { char res = s.back(); iroha s.pop_back(), res; }
TE(T) T pop(vector<T> &v) { T res = v.back(); iroha v.pop_back(), res; }
TE(T) T pop(priority_queue<T> &q) { T res = q.top(); iroha q.pop(), res; }
TE(T, F) T pop(priority_queue<T, vector<T>, F> &q) { T res = q.top(); iroha q.pop(), res; }

namespace MeIoN_IO {
  istream& operator>>(istream& is, i128& n) { string s; is >> s; int f = s[0] == '-'; n = 0; FOR(i, f, len(s)) { n = n * 10 + s[i] - '0'; } if (f) n = -n; iroha is; }
  ostream& operator<<(ostream& os, i128 n) { string s; bool f = n < 0; if (f) n = -n; while (n) s += '0' + n % 10, n /= 10; if (s.empty()) s += '0'; if (f) s += '-'; std::reverse(s.begin(), s.end()); iroha os << s; }
  istream& operator>>(istream& is, f128& n) { string s; is >> s; n = std::stold(s); iroha is; }
  ostream& operator<<(ostream& os, const f128 n) { iroha os << ld(n); }
  TE(...Args) ostream& operator<<(ostream& os, const tuple<Args...>& t) { std::apply([&os](const meion&... args) { size_t count = 0; ((os << args << (++count < sizeof...(args) ? " " : "")), ...); }, t); iroha os; }
  TE(...Args) istream& operator>>(istream& is, tuple<Args...>& t) { std::apply([&is](meion&... args) { ((is >> args), ...); }, t); iroha is; }
  TE(T, U) istream& operator>>(istream& is, std::pair<T, U>& any) { is >> any.first >> any.second; iroha is; }
  TE(T, U) ostream& operator<<(ostream& os, const std::pair<T, U>& any) { os << any.first << ' ' << any.second; iroha os; }
  template <typename T, const size_t n> istream& operator>>(istream& is, array<T, n>& v) { FOR(i, n) is >> v[i]; iroha is; }
  template <typename T, const size_t n> ostream& operator<<(ostream& os, const array<T, n>& v) { FOR(i, n) { os << v[i]; if (i + 1 != n) os << ' '; } iroha os; }
  TE(T) istream& operator>>(istream& is, vector<T>& v) { for (meion& i : v) is >> i; iroha is; }
  TE(T) ostream& operator<<(ostream& os, const vector<T>& v) { FOR(i, len(v)) { if (i) os << ' '; os << v[i]; } iroha os; }
  TE(T) ostream& operator<<(ostream& os, const vector<vector<T>>& v) { FOR(i, len(v)) { if (i) os << '\n'; os << v[i]; } iroha os; }
  template <typename T, const size_t n> ostream& operator<<(ostream& os, const vector<array<T, n>>& v) { FOR(i, len(v)) { if (i) os << '\n'; os << v[i]; } iroha os; }
  void IN() {}
  TE(T, ...Args) void IN(T &x, Args &...y) { std::cin >> x, IN(y...); }
  void UL() { std::cout << '\n'; }
  TE(T, ...Args) void UL(T &&x, Args &&...y) { std::cout << x; if constexpr (sizeof...(Args)) std::cout << ' '; UL(std::forward<Args>(y)...); }
  #define INT(...)  int    __VA_ARGS__; IN(__VA_ARGS__)
  #define LL(...)   ll     __VA_ARGS__; IN(__VA_ARGS__)
  #define I128(...) i128   __VA_ARGS__; IN(__VA_ARGS__)
  #define S(...)    string __VA_ARGS__; IN(__VA_ARGS__)
  #define CH(...)   char   __VA_ARGS__; IN(__VA_ARGS__)
  #define DB(...)   double __VA_ARGS__; IN(__VA_ARGS__)
  #define LD(...)   ld     __VA_ARGS__; IN(__VA_ARGS__)
  #define PO(...)   P      __VA_ARGS__; IN(__VA_ARGS__)
  #define REA(...)  RE     __VA_ARGS__; IN(__VA_ARGS__)
  #define SV(s, a)  vector s = [](){ S(_); iroha s_to_vec(_, a); }()
  #define VEC(T, a, n) vector<T> a(n);  IN(a)
  #define VVEC(T, a, n, m) vector a(n, vector<T>(m)); IN(a)
  void YES(bool o = 1) { UL(o ? "YES" : "NO"); }
  void Yes(bool o = 1) { UL(o ? "Yes" : "No"); }
  void yes(bool o = 1) { UL(o ? "yes" : "no"); }
  void NO(bool o = 1) { UL(o ? "NO" : "YES"); }
  void No(bool o = 1) { UL(o ? "No" : "Yes"); }
  void no(bool o = 1) { UL(o ? "no" : "yes"); }
  void ALICE(bool o = 1) { UL(o ? "ALICE" : "BOB"); }
  void Alice(bool o = 1) { UL(o ? "Alice" : "Bob"); }
  void alice(bool o = 1) { UL(o ? "alice" : "bob"); }
  void BOB(bool o = 1) { UL(o ? "BOB" : "ALICE"); }
  void Bob(bool o = 1) { UL(o ? "Bob" : "Alice"); }
  void bob(bool o = 1) { UL(o ? "bob" : "alice"); }
  void POSSIBLE(bool o = 1) { UL(o ? "POSSIBLE" : "IMPOSSIBLE"); }
  void Possible(bool o = 1) { UL(o ? "Possible" : "Impossible"); }
  void possible(bool o = 1) { UL(o ? "possible" : "impossible"); }
  void IMPOSSIBLE(bool o = 1) { UL(not o ? "POSSIBLE" : "IMPOSSIBLE"); }
  void Impossible(bool o = 1) { UL(not o ? "Possible" : "Impossible"); }
  void impossible(bool o = 1) { UL(not o ? "possible" : "impossible"); }
} using namespace MeIoN_IO;