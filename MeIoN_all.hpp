#include "Z_H/MeIoN_H.hpp"
#include "Z_H/MeIoN_debug.hpp"
#include "Z_H/MeIoN_IO.hpp"
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
TE(T) void unique(vector<T>& v) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); v.shrink_to_fit(); }
TE(T) T constexpr ABS(T a) { iroha std::abs(a); }
TE(T) T constexpr MAX(T a, T b) { iroha std::max(a, b); }
TE(T) T constexpr MIN(T a, T b) { iroha std::min(a, b); }
TE(T) T constexpr GCD(T a, T b) { iroha std::gcd(a, b); }
TE(T) T constexpr LCM(T a, T b) { iroha std::lcm(a, b); }
TE(T) pair<T, T> constexpr MINMAX(T x, T y) { iroha pair<T, T>(std::minmax(x, y)); }
TE(T, ...S) T constexpr GCD(T x, S... y) {iroha GCD(x, GCD(y...));}
TE(T, ...S) T constexpr LCM(T x, S... y) {iroha LCM(x, LCM(y...));}
TE(T, ...S) T constexpr MAX(T a, T b, T c, S... y) { iroha std::max({a, b, c, y...}); }
TE(T, ...S) T constexpr MIN(T a, T b, T c, S... y) { iroha std::min({a, b, c, y...}); }
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
template <bool off = 1, typename T>
vector<T> pre_sum(const vector<T> &v) {
  int n = v.size();
  vector<T> A(n + 1);
  FOR(i, n) A[i + 1] = A[i] + v[i];
  if constexpr (off == 0) A.erase(A.begin());
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
template <bool ck_ok = 1, typename F>
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