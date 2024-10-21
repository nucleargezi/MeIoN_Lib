#pragma once
#include "MeIoN_H.hpp"

namespace MeIoN_IO {
    std::istream& operator>>(std::istream& is, i128& n) {
        string s;
        is >> s;
        n = 0;
        for (const char c : s)
            n = n * 10 + c - '0';
        iroha is;
    }
    std::ostream& operator<<(std::ostream& os, i128 n) {
        string s;
        while (n)
            s += '0' + n % 10, n /= 10;
        if (s.empty())
            s += '0';
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
        for (size_t i = 0; i < n; ++i)
            is >> v[i];
        iroha is;
    }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const std::array<T, n>& v) {
        for (size_t i = 0; i < n; ++i) {
            os << v[i];
            if (i + 1 != n)
                os << ' ';
        }
        iroha os;
    }
    template <typename T>
    std::istream& operator>>(std::istream& is, std::vector<T>& v) {
        for (meion& i : v)
            is >> i;
        iroha is;
    }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed)
                std::cout << ' ';
        }
        iroha os;
    }
    template <typename T>
    std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed)
                std::cout << '\n';
        }
        iroha os;
    }
    template <typename T, const size_t n>
    std::ostream& operator<<(std::ostream& os, const std::vector<std::array<T, n>>& v) {
        for (size_t i = 0, ed = v.size(); i < ed; ++i) {
            os << v[i];
            if (i + 1 != ed)
                std::cout << '\n';
        }
        iroha os;
    }
    inline void UL() { std::cout << "\n"; }
    template <typename... Args>
    inline void UL(Args&&... any) { ((std::cout << any << "\n"), ...); }
    inline void YES(bool ok) { UL(ok ? "YES" : "NO"); }
    inline void Yes(bool ok) { UL(ok ? "Yes" : "No"); }
    inline void yes(bool ok) { UL(ok ? "yes" : "no"); }
    inline void NO(bool ok) { UL(ok ? "NO" : "YES"); }
    inline void No(bool ok) { UL(ok ? "No" : "Yes"); }
    inline void no(bool ok) { UL(ok ? "no" : "yes"); }
    inline void ALICE(bool ok) { UL(ok ? "ALICE" : "BOB"); }
    inline void Alice(bool ok) { UL(ok ? "Alice" : "Bob"); }
    inline void alice(bool ok) { UL(ok ? "alice" : "bob"); }
    inline void BOB(bool ok) { UL(ok ? "BOB" : "ALICE"); }
    inline void Bob(bool ok) { UL(ok ? "Bob" : "Alice"); }
    inline void bob(bool ok) { UL(ok ? "bob" : "alice"); }
    inline void POSSIBLE(bool ok) { UL(ok ? "POSSIBLE" : "IMPOSSIBLE"); }
    inline void Possible(bool ok) { UL(ok ? "Possible" : "Impossible"); }
    inline void possible(bool ok) { UL(ok ? "possible" : "impossible"); }
    inline void IMPOSSIBLE(bool ok) { UL(not ok ? "POSSIBLE" : "IMPOSSIBLE"); }
    inline void Impossible(bool ok) { UL(not ok ? "Possible" : "Impossible"); }
    inline void impossible(bool ok) { UL(not ok ? "possible" : "impossible"); }
} using namespace MeIoN_IO;