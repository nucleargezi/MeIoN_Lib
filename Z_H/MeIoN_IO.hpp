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
    std::istream& operator>>(std::istream& is, std::array<T, n>& v) {
        for (size_t i = 0; i < n; ++i) is >> v[i];
        iroha is;
    }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const std::array<T, n>& v) {
        for (size_t i = 0; i < n; ++i) {
            os << v[i];
            if (i + 1 != n) os << ' ';
        }
        iroha os;
    }
    template <typename T>
    std::istream& operator>>(std::istream& is, std::vector<T>& v) {
        for (meion& i : v) is >> i;
        iroha is;
    }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed) std::cout << ' ';
        }
        iroha os;
    }
    template <typename T>
    std::ostream& operator<<(std::ostream& os,
                            const std::vector<std::vector<T>>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed) std::cout << '\n';
        }
        iroha os;
    }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os,
                            const std::vector<std::array<T, n>>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed) std::cout << '\n';
        }
        iroha os;
    }
    inline void YES(bool ok = true) { std::cout << (ok ? "YES" : "NO") << '\n'; }
    inline void Yes(bool ok = true) { std::cout << (ok ? "Yes" : "No") << '\n'; }
    inline void yes(bool ok = true) { std::cout << (ok ? "yes" : "no") << '\n'; }
    inline void NO(bool ok = true) { std::cout << (ok ? "NO" : "YES") << '\n'; }
    inline void No(bool ok = true) { std::cout << (ok ? "No" : "Yes") << '\n'; }
    inline void no(bool ok = true) { std::cout << (ok ? "no" : "yes") << '\n'; }
    inline void ALICE(bool ok = true) { std::cout << (ok ? "ALICE" : "BOB") << '\n'; }
    inline void Alice(bool ok = true) { std::cout << (ok ? "Alice" : "Bob") << '\n'; }
    inline void alice(bool ok = true) { std::cout << (ok ? "alice" : "bob") << '\n'; }
    inline void BOB(bool ok = true) { std::cout << (ok ? "BOB" : "ALICE") << '\n'; }
    inline void Bob(bool ok = true) { std::cout << (ok ? "Bob" : "Alice") << '\n'; }
    inline void bob(bool ok = true) { std::cout << (ok ? "bob" : "alice") << '\n'; }
    inline void POSSIBLE(bool ok = true) { std::cout << (ok ? "POSSIBLE" : "IMPOSSIBLE") << '\n'; }
    inline void Possible(bool ok = true) { std::cout << (ok ? "Possible" : "Impossible") << '\n'; }
    inline void possible(bool ok = true) { std::cout << (ok ? "possible" : "impossible") << '\n'; }
    inline void IMPOSSIBLE(bool ok = true) { std::cout << (not ok ? "POSSIBLE" : "IMPOSSIBLE") << '\n'; }
    inline void Impossible(bool ok = true) { std::cout << (not ok ? "Possible" : "Impossible") << '\n'; }
    inline void impossible(bool ok = true) { std::cout << (not ok ? "possible" : "impossible") << '\n'; }
} using namespace MeIoN_IO;