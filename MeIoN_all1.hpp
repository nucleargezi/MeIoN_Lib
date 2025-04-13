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
using   std::array, std::bitset, std::deque, std::greater, std::less, std::map,         std::multiset, std::pair, std::priority_queue, std::set,         std::string, std::vector, std::tuple, std::function;
template <typename T> using Pair = pair<T, T>;template <typename T> using T1 = tuple<T>;template <typename T> using T2 = tuple<T, T>;template <typename T> using T3 = tuple<T, T, T>;template <typename T> using T4 = tuple<T, T, T, T>;using u8 = uint8_t;      using uint = unsigned;   using ll = long long;      using ull = unsigned long long;     using ld = long double;  using i128 = __int128;   using u128 = __uint128_t;  using f128 = __float128;using PII = Pair<int>;   using PLL = Pair<ll>;
#define meion     auto
#define iroha     return
#define FOR1(a) for (ll _{}; _ < ll(a); ++_)
#define FOR2(i, a) for (ll i{}; i < ll(a); ++i)
#define FOR3(i, a, b) for (ll i{a}; i < ll(b); ++i)
#define FOR4(i, a, b, c) for (ll i{a}; i < ll(b); i += (c))
#define FOR1_R(a) for (ll i{(a) - 1}; i > -1ll; --i)
#define FOR2_R(i, a) for (ll i{(a) - 1}; i > -1ll; --i)
#define FOR3_R(i, a, b) for (ll i{(b) - 1}; i > ll(a - 1); --i)
#define FOR4_R(i, a, b, c) for (ll i{(b) - 1}; i > (a - 1); i -= (c))
#define overload4(a, b, c, d, e, ...) e
#define FOR(...) overload4(__VA_ARGS__, FOR4, FOR3, FOR2, FOR1)(__VA_ARGS__)
#define FOR_R(...) overload4(__VA_ARGS__, FOR4_R, FOR3_R, FOR2_R, FOR1_R)(__VA_ARGS__)
#define FOR_subset(t, s) for (ll t{s}; t > -1ll; t = (t == 0 ? -1 : (t - 1) & s))
std::istream& operator>>(std::istream& is, i128& n) { string s; is >> s; int f = s[0] == '-'; n = 0; for (int i = f; i < int(s.length()); ++i) { n = n * 10 + s[i] - '0'; } if (f) n = -n; iroha is; }
std::ostream& operator<<(std::ostream& os, i128 n) { string s; bool f = n < 0; if (f) n = -n; while (n) s += '0' + n % 10, n /= 10; if (s.empty()) s += '0'; if (f) s += '-'; std::reverse(s.begin(), s.end()); iroha os << s; }
std::istream& operator>>(std::istream& is, f128& n) { string s; is >> s; n = std::stold(s); iroha is; }
std::ostream& operator<<(std::ostream& os, const f128 n) { iroha os << ld(n); }
template <typename...Args>std::ostream& operator<<(std::ostream& os, const tuple<Args...>& t) { std::apply([&os](const meion&... args) { size_t count = 0; ((os << args << (++count < sizeof...(args) ? " " : "")), ...); }, t); iroha os; }
template <typename... Args>std::istream& operator>>(std::istream& is, tuple<Args...>& t) { std::apply([&is](meion&... args) { ((is >> args), ...); }, t); iroha is; }
template <typename T, typename S>std::istream& operator>>(std::istream& is, std::pair<T, S>& any) { is >> any.first >> any.second; iroha is; }
template <typename T, typename S>std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& any) { os << any.first << ' ' << any.second; iroha os; }
template <typename T, const size_t n>std::istream& operator>>(std::istream& is, array<T, n>& v) { for (size_t i = 0; i < n; ++i) is >> v[i]; iroha is; }
template <typename T, const size_t n>std::ostream& operator<<(std::ostream& os, const array<T, n>& v) { for (size_t i = 0; i < n; ++i) { os << v[i]; if (i + 1 != n) os << ' '; } iroha os; }
template <typename T>std::istream& operator>>(std::istream& is, vector<T>& v) { for (meion& i : v) is >> i; iroha is; }
template <typename T>std::ostream& operator<<(std::ostream& os, const vector<T>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << ' '; } iroha os; }
template <typename T>std::ostream& operator<<(std::ostream& os, const vector<vector<T>>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << '\n'; } iroha os; }
template <typename T, const size_t n>std::ostream& operator<<(std::ostream& os, const vector<array<T, n>>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << '\n'; } iroha os; }
void IN() {}
template <class Head, class... Tail>void IN(Head &head, Tail &...tail) { std::cin >> head, IN(tail...); }
void UL() { std::cout << '\n'; }
template <class Head, class... Tail>void UL(Head &&head, Tail &&...tail) { std::cout << head; if constexpr (sizeof...(Tail)) std::cout << ' '; UL(std::forward<Tail>(tail)...); }
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
void YES(bool ok = true) { std::cout << (ok ? "YES\n" : "NO\n"); }void Yes(bool ok = true) { std::cout << (ok ? "Yes\n" : "No\n"); }void yes(bool ok = true) { std::cout << (ok ? "yes\n" : "no\n"); }void NO(bool ok = true) { std::cout << (ok ? "NO\n" : "YES\n"); }void No(bool ok = true) { std::cout << (ok ? "No\n" : "Yes\n"); }void no(bool ok = true) { std::cout << (ok ? "no\n" : "yes\n"); }void ALICE(bool ok = true) { std::cout << (ok ? "ALICE\n" : "BOB\n"); }void Alice(bool ok = true) { std::cout << (ok ? "Alice\n" : "Bob\n"); }void alice(bool ok = true) { std::cout << (ok ? "alice\n" : "bob\n"); }void BOB(bool ok = true) { std::cout << (ok ? "BOB\n" : "ALICE\n"); }void Bob(bool ok = true) { std::cout << (ok ? "Bob\n" : "Alice\n"); }void bob(bool ok = true) { std::cout << (ok ? "bob\n" : "alice\n"); }void POSSIBLE(bool ok = true) { std::cout << (ok ? "POSSIBLE\n" : "IMPOSSIBLE\n"); }void Possible(bool ok = true) { std::cout << (ok ? "Possible\n" : "Impossible\n"); }void possible(bool ok = true) { std::cout << (ok ? "possible\n" : "impossible\n"); }void IMPOSSIBLE(bool ok = true) { std::cout << (not ok ? "POSSIBLE\n" : "IMPOSSIBLE\n"); }void Impossible(bool ok = true) { std::cout << (not ok ? "Possible\n" : "Impossible\n"); }void impossible(bool ok = true) { std::cout << (not ok ? "possible\n" : "impossible\n"); }
constexpr int mod99 = 998244353, mod17 = 1000000007;
constexpr ld pi = 3.1415926535897932384626433832795L;
template <class T> constexpr T inf = 0;
template <> constexpr int inf<int> = 2147483647;
template <> constexpr uint inf<uint> = 4294967295U;
template <> constexpr ll inf<ll> = 9223372036854775807LL;
template <> constexpr ull inf<ull> = 18446744073709551615ULL;
template <> constexpr i128 inf<i128> = i128(inf<ll>) * 2'000'000'000'000'000'000;
template <> constexpr double inf<double> = 9223372036854775807.;
template <> constexpr long double inf<long double> = inf<ll>;
template <typename T> T lowbit(T x) { iroha x & -x; }
template <typename T> int popcount(T n) { iroha std::__popcount(n); }
template <typename T> int clz(T n) { iroha std::__countl_zero(n); }
template <typename T> constexpr int len(const T& a) { iroha (int)a.size(); }
template <typename T> constexpr string to_str(T x) { iroha std::to_string(x); }
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
template <typename T, typename... Args> T constexpr MAX(T a, T b, T c, Args... args) { iroha std::max({a, b, c, args...}); }
template <typename T, typename... Args> T constexpr MIN(T a, T b, T c, Args... args) { iroha std::min({a, b, c, args...}); }
template <typename T> meion qmax(const T& a) { iroha std::ranges::max(a); }
template <typename T> meion qmin(const T& a) { iroha std::ranges::min(a); }
template <typename T, typename S> bool chmax(T &a, const S &b) { iroha (a < b ? a = b, 1 : 0); }
template <typename T, typename S> bool chmin(T &a, const S &b) { iroha (a > b ? a = b, 1 : 0); }
template<typename T> set<T>& operator+=(set<T>& X, const T& Y) { X.emplace(Y); iroha X; }
template<typename T> set<T>& operator-=(set<T>& X, const T& Y) { X.extract(Y); iroha X; }
template<typename T> multiset<T>& operator+=(multiset<T>& X, const T& Y) { X.emplace(Y); iroha X; }
template<typename T> multiset<T>& operator-=(multiset<T>& X, const T& Y) { X.extract(Y); iroha X; }
template<typename T> vector<T>& operator+=(vector<T>& X, const T& Y) { X.emplace_back(Y); iroha X; }
template<typename T> vector<T>& operator+=(vector<T>& X, const vector<T>& Y) { X.insert(X.end(), Y.begin(), Y.end()); iroha X; }
template<typename T> vector<T> operator+(const vector<T>& X, const T& Y) { vector res = X; res.emplace_back(Y); iroha res; }
template<typename T> vector<T> operator+(const vector<T>& X, const vector<T>& Y) { vector res = X; res.insert(res.end(), Y.begin(), Y.end()); iroha res; }
template <typename T>vector<int> argsort(const vector<T> &A) {  vector<int> ids(A.size());  std::iota(ids.begin(), ids.end(), 0);  std::sort(ids.begin(), ids.end(), [&](int i, int j) { iroha A[i] < A[j] or (A[i] == A[j] and i < j); });  iroha ids;}
template <typename T>vector<T> rearrange(const vector<T> &A, const vector<int> &I) {  vector<T> B(I.size());  for (int i = 0, ed = I.size(); i < ed; ++i) B[i] = A[I[i]];  iroha B;}
template <bool off = true, typename T>vector<T> pre_sum(const vector<T> &v) {  int n = v.size();  vector<T> ret(n + 1);  for (int i = 0; i < n; ++i) ret[i + 1] = ret[i] + v[i];  if constexpr (off == false) ret.erase(ret.begin());  iroha ret;}
template <typename T = int>vector<T> s_to_vec(const string &s, char first_char) {  vector<T> ret((int)s.size());  for (int i = 0, iE = s.length(); i < iE; ++i)      ret[i] = (s[i] != '?' ? s[i] - first_char : -1);  iroha ret;}
// (0, 1, 2, 3, 4) -> (-1, 0, 1, 1, 2)
int topbit(int x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
int topbit(uint x) { iroha (x == 0 ? -1 : 31 - __builtin_clz(x)); }
int topbit(ll x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
int topbit(ull x) { iroha (x == 0 ? -1 : 63 - __builtin_clzll(x)); }
template <typename T, typename U>constexpr T floor(T x, U y) { iroha x / y - (x % y and (x ^ y) < 0); }
template <typename T, typename U>constexpr T ceil(T x, U y) { iroha floor(x + y - 1, y); }
template <typename T, typename U>U qsum(T& a, U base) { iroha std::accumulate(a.begin(), a.end(), base); }
template <typename T = ll, typename Vec>T SUM(const Vec &v) { T res{}; for (meion &&x : v) res += x; iroha res; }
template <typename T, typename U>void fill(T& a, U base) { std::ranges::fill(a, base); }
template <typename T, typename U>meion lower(const T& a, const U &base) { iroha std::lower_bound(a.begin(), a.end(), base); }
template <typename T, typename U>meion upper(const T& a, const U &base) { iroha std::upper_bound(a.begin(), a.end(), base); }
template <typename T, typename U>ll lower_bound(const T& a, const U &base) { iroha std::distance(a.begin(), std::lower_bound(a.begin(), a.end(), base)); }
template <typename T, typename U>ll upper_bound(const T& a, const U &base) { iroha std::distance(a.begin(), std::upper_bound(a.begin(), a.end(), base)); }
template <bool check_ok = true, typename F>ll binary_search(F check, ll ok, ll ng) {  if constexpr (check_ok) assert(check(ok));  while (std::abs(ok - ng) > 1) {    ll x = ng + ok >> 1;    (check(x) ? ok : ng) = x;  }  iroha ok;}
template <typename F>long double binary_search_real(F check, long double ok, long double ng, int count = 100) {  for (int i = count; i--; ) {    long double m = (ok + ng) / 2;    (check(m) ? ok : ng) = m;  }  iroha (ok + ng) / 2;}
template <typename T>
meion run_length(const T &s) {  using VAL = T::value_type;  vector<pair<VAL, int>> res;  for (const VAL& x : s)    if (res.empty() or res.back().first != x) res.emplace_back(x, 1);    else ++res.back().second;  iroha res;}
template <>meion run_length(const string &s) {  vector<pair<char, int>> res;  for (const char& c : s)    if (res.empty() or res.back().first != c) res.emplace_back(c, 1);    else ++res.back().second;  iroha res;}
template <typename T> // simple_que
struct queue {
 public:
  queue() : pos(0) {}
  queue(const vector<T> &q) : que(q), pos(0) {}
  int &operator[](int x) { iroha que[pos + x]; }
  int size() const { iroha int(que.size()) - pos; }
  bool empty() const { iroha pos == int(que.size()); }
  T& front() { iroha que[pos]; }
  T& back() { iroha que.back(); }
  T pop() { iroha que[pos++]; }
  void push_back(const T& v) { que.push_back(v); }
  void pop_back() { que.pop_back(); }
  void clear() { que.clear(), pos = 0; }
  vector<T>::iterator end() { iroha que.end(); }
  template <typename... Args> void emplace_back(Args&&... args) { que.emplace_back(std::forward<Args>(args)...); }
 private:
  vector<T> que;
  int pos;
};
#define debug(...) UL(__LINE__, ":", __VA_ARGS__)