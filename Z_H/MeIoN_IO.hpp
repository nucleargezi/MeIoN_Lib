istream& operator>>(istream& I, i128& n) { string s; I >> s; int f = s[0] == '-'; n = 0; FOR(i, f, s.size()) { n = n * 10 + s[i] - '0'; } if (f) n = -n; iroha I; }
ostream& operator<<(ostream& O, i128 n) { string s; bool f = n < 0; if (f) n = -n; while (n) s += '0' + n % 10, n /= 10; if (s.empty()) s += '0'; if (f) s += '-'; std::reverse(s.begin(), s.end()); iroha O << s; }
istream& operator>>(istream& I, f128& n) { string s; I >> s; n = std::stold(s); iroha I; }
ostream& operator<<(ostream& O, const f128 n) { iroha O << ld(n); }
TE(...S) istream& operator>>(istream& I, tuple<S...>& t) { std::apply([&I](meion&... args) { ((I >> args), ...); }, t); iroha I; }
TE(T, U) istream& operator>>(istream& I, pair<T, U>& x) { I >> x.first >> x.second; iroha I; }
TE(T, U) ostream& operator<<(ostream& O, const pair<T, U>& x) { O << x.first << ' ' << x.second; iroha O; }
template <typename T, const size_t n> istream& operator>>(istream& I, array<T, n>& v) { FOR(i, n) I >> v[i]; iroha I; }
template <typename T, const size_t n> ostream& operator<<(ostream& O, const array<T, n>& v) { FOR(i, n) { O << v[i]; if (i + 1 != n) O << ' '; } iroha O; }
TE(T) istream& operator>>(istream& I, vector<T>& v) { for (meion& i : v) I >> i; iroha I; }
TE(T) ostream& operator<<(ostream& O, const vector<T>& v) { FOR(i, v.size()) { if (i) O << ' '; O << v[i]; } iroha O; }
TE(T) ostream& operator<<(ostream& O, const vector<vector<T>>& v) { FOR(i, v.size()) { if (i) O << '\n'; O << v[i]; } iroha O; }
template <typename T, const size_t n> ostream& operator<<(ostream& O, const vector<array<T, n>>& v) { FOR(i, v.size()) { if (i) O << '\n'; O << v[i]; } iroha O; }
void IN() {}
TE(T, ...S) void IN(T &x, S &...y) { std::cin >> x, IN(y...); }
void UL() { std::cout << '\n'; }
TE(T, ...S) void UL(T &&x, S &&...y) { std::cout << x; if constexpr (sizeof...(S)) std::cout << ' '; UL(std::forward<S>(y)...); }
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
void YES(bool o = 1) { UL(o ? "YES" : "NO"); } void NO(bool o = 1) { YES(not o); }
void Yes(bool o = 1) { UL(o ? "Yes" : "No"); } void No(bool o = 1) { Yes(not o); }
void yes(bool o = 1) { UL(o ? "yes" : "no"); } void no(bool o = 1) { yes(not o); }
void ALICE(bool o = 1) { UL(o ? "ALICE" : "BOB"); } void BOB(bool o = 1) { ALICE(not o); }
void Alice(bool o = 1) { UL(o ? "Alice" : "Bob"); } void Bob(bool o = 1) { Alice(not o); }
void alice(bool o = 1) { UL(o ? "alice" : "bob"); } void bob(bool o = 1) { alice(not o); }
void POSSIBLE(bool o = 1) { UL(o ? "POSSIBLE" : "IMPOSSIBLE"); } void IMPOSSIBLE(bool o = 1) { POSSIBLE(not o); }
void Possible(bool o = 1) { UL(o ? "Possible" : "Impossible"); } void Impossible(bool o = 1) { Possible(not o); }
void possible(bool o = 1) { UL(o ? "possible" : "impossible"); } void impossible(bool o = 1) { possible(not o); }