namespace MeIoN_IO {
  istream& operator>>(istream& is, i128& n) { string s; is >> s; int f = s[0] == '-'; n = 0; for (int i = f; i < int(s.length()); ++i) { n = n * 10 + s[i] - '0'; } if (f) n = -n; iroha is; }
  ostream& operator<<(ostream& os, i128 n) { string s; bool f = n < 0; if (f) n = -n; while (n) s += '0' + n % 10, n /= 10; if (s.empty()) s += '0'; if (f) s += '-'; std::reverse(s.begin(), s.end()); iroha os << s; }
  istream& operator>>(istream& is, f128& n) { string s; is >> s; n = std::stold(s); iroha is; }
  ostream& operator<<(ostream& os, const f128 n) { iroha os << ld(n); }
  TE(...Args) ostream& operator<<(ostream& os, const tuple<Args...>& t) { std::apply([&os](const meion&... args) { size_t count = 0; ((os << args << (++count < sizeof...(args) ? " " : "")), ...); }, t); iroha os; }
  TE(...Args) istream& operator>>(istream& is, tuple<Args...>& t) { std::apply([&is](meion&... args) { ((is >> args), ...); }, t); iroha is; }
  TE(T, U) istream& operator>>(istream& is, std::pair<T, U>& any) { is >> any.first >> any.second; iroha is; }
  TE(T, U) ostream& operator<<(ostream& os, const std::pair<T, U>& any) { os << any.first << ' ' << any.second; iroha os; }
  template <typename T, const size_t n> istream& operator>>(istream& is, array<T, n>& v) { for (size_t i = 0; i < n; ++i) is >> v[i]; iroha is; }
  template <typename T, const size_t n> ostream& operator<<(ostream& os, const array<T, n>& v) { for (size_t i = 0; i < n; ++i) { os << v[i]; if (i + 1 != n) os << ' '; } iroha os; }
  TE(T) istream& operator>>(istream& is, vector<T>& v) { for (meion& i : v) is >> i; iroha is; }
  TE(T) ostream& operator<<(ostream& os, const vector<T>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << ' '; } iroha os; }
  TE(T) ostream& operator<<(ostream& os, const vector<vector<T>>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << '\n'; } iroha os; }
  template <typename T, const size_t n> ostream& operator<<(ostream& os, const vector<array<T, n>>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << '\n'; } iroha os; }
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
  void YES(bool ok = 1) { UL(ok ? "YES" : "NO"); }
  void Yes(bool ok = 1) { UL(ok ? "Yes" : "No"); }
  void yes(bool ok = 1) { UL(ok ? "yes" : "no"); }
  void NO(bool ok = 1) { UL(ok ? "NO" : "YES"); }
  void No(bool ok = 1) { UL(ok ? "No" : "Yes"); }
  void no(bool ok = 1) { UL(ok ? "no" : "yes"); }
  void ALICE(bool ok = 1) { UL(ok ? "ALICE" : "BOB"); }
  void Alice(bool ok = 1) { UL(ok ? "Alice" : "Bob"); }
  void alice(bool ok = 1) { UL(ok ? "alice" : "bob"); }
  void BOB(bool ok = 1) { UL(ok ? "BOB" : "ALICE"); }
  void Bob(bool ok = 1) { UL(ok ? "Bob" : "Alice"); }
  void bob(bool ok = 1) { UL(ok ? "bob" : "alice"); }
  void POSSIBLE(bool ok = 1) { UL(ok ? "POSSIBLE" : "IMPOSSIBLE"); }
  void Possible(bool ok = 1) { UL(ok ? "Possible" : "Impossible"); }
  void possible(bool ok = 1) { UL(ok ? "possible" : "impossible"); }
  void IMPOSSIBLE(bool ok = 1) { UL(not ok ? "POSSIBLE" : "IMPOSSIBLE"); }
  void Impossible(bool ok = 1) { UL(not ok ? "Possible" : "Impossible"); }
  void impossible(bool ok = 1) { UL(not ok ? "possible" : "impossible"); }
} using namespace MeIoN_IO;