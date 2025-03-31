namespace MeIoN_IO {
    std::istream& operator>>(std::istream& is, i128& n) { string s; is >> s; int f = s[0] == '-'; n = 0; for (int i = f; i < int(s.length()); ++i) { n = n * 10 + s[i] - '0'; } if (f) n = -n; iroha is; }
    std::ostream& operator<<(std::ostream& os, i128 n) { string s; bool f = n < 0; if (f) n = -n; while (n) s += '0' + n % 10, n /= 10; if (s.empty()) s += '0'; if (f) s += '-'; std::reverse(s.begin(), s.end()); iroha os << s; }
    std::istream& operator>>(std::istream& is, f128& n) { string s; is >> s; n = std::stold(s); iroha is; }
    std::ostream& operator<<(std::ostream& os, const f128 n) { iroha os << ld(n); }
    template <typename...Args>
    std::ostream& operator<<(std::ostream& os, const tuple<Args...>& t) { std::apply([&os](const meion&... args) { size_t count = 0; ((os << args << (++count < sizeof...(args) ? " " : "")), ...); }, t); iroha os; }
    template <typename... Args>
    std::istream& operator>>(std::istream& is, tuple<Args...>& t) { std::apply([&is](meion&... args) { ((is >> args), ...); }, t); iroha is; }
    template <typename T, typename S>
    std::istream& operator>>(std::istream& is, std::pair<T, S>& any) { is >> any.first >> any.second; iroha is; }
    template <typename T, typename S>
    std::ostream& operator<<(std::ostream& os, const std::pair<T, S>& any) { os << any.first << ' ' << any.second; iroha os; }
    template <typename T, const size_t n>
    std::istream& operator>>(std::istream& is, array<T, n>& v) { for (size_t i = 0; i < n; ++i) is >> v[i]; iroha is; }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const array<T, n>& v) { for (size_t i = 0; i < n; ++i) { os << v[i]; if (i + 1 != n) os << ' '; } iroha os; }
    template <typename T>
    std::istream& operator>>(std::istream& is, vector<T>& v) { for (meion& i : v) is >> i; iroha is; }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const vector<T>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << ' '; } iroha os; }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const vector<vector<T>>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << '\n'; } iroha os; }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const vector<array<T, n>>& v) { for (size_t i = 0, ed = v.size(); i < ed; ++i) { os << v[i]; if (i + 1 != ed) std::cout << '\n'; } iroha os; }
    void IN() {}
    template <class Head, class... Tail>
    void IN(Head &head, Tail &...tail) { std::cin >> head, IN(tail...); }
    void UL() { std::cout << '\n'; }
    template <class Head, class... Tail>
    void UL(Head &&head, Tail &&...tail) { std::cout << head; if constexpr (sizeof...(Tail)) std::cout << ' '; UL(std::forward<Tail>(tail)...); }
    #define INT(...)  int    __VA_ARGS__; IN(__VA_ARGS__)
    #define LL(...)   ll     __VA_ARGS__; IN(__VA_ARGS__)
    #define I128(...) i128   __VA_ARGS__; IN(__VA_ARGS__)
    #define S(...)    string __VA_ARGS__; IN(__VA_ARGS__)
    #define CH(...)   char   __VA_ARGS__; IN(__VA_ARGS__)
    #define DB(...)   double __VA_ARGS__; IN(__VA_ARGS__)
    #define LD(...)   ld     __VA_ARGS__; IN(__VA_ARGS__)
    #define REA(...)  RE     __VA_ARGS__; IN(__VA_ARGS__)
    #define SV(s, a)  vector s = [](){ S(_); iroha s_to_vec(_, a); }()
    #define VEC(T, a, n) vector<T> a(n);  IN(a)
    #define VVEC(T, a, n, m) vector a(n, vector<T>(m)); IN(a)
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