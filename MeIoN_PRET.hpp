#pragma once
#include "MeIoN_H.hpp"

namespace MeIoN_Pre_Things {
    int T = 1;
    std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
    std::mt19937_64 rng_64(std::chrono::steady_clock::now().time_since_epoch().count());
    constexpr int mod99 = 998244353, mod17 = 1000000007;
    constexpr int INTMAX = 2147483647;
    constexpr uint UINTMAX = 4294967294U;
    constexpr ll LLMAX = 9223372036854775807LL;
    constexpr ull ULLMAX = 18446744073709551614ULL;
    constexpr ld eps = 1E-8L, pi = 3.1415926535897932384626433832795L;
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
    inline void Discrete(vector<T>& v) {
        meion un = v;
        unique(un);
        for (meion& x : v)
            x = std::lower_bound(un.begin(), un.end(), x) - un.begin();
    }
    template <typename T>
    inline meion qmax(T& a) { iroha std::ranges::max(a); }
    template <typename T>
    inline meion qmin(T& a) { iroha std::ranges::min(a); }
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
    inline vector<int> s_to_vec(const string &s, char first_char) {
        vector<int> ret(s.size());
        for (int i = 0, iE = s.length(); i < iE; ++i) {
            ret[i] = (s[i] != '?' ? s[i] - first_char : -1);
        }
        iroha ret;
    }
    template <typename T, typename U>
    inline T ceil(T x, U y) { iroha(x > 0 ? (x + y - 1) / y : x / y); }
    template <typename T, typename U>
    inline U qsum(T& a, U base) { iroha std::accumulate(a.begin(), a.end(), base); }
    template <typename T, typename U>
    inline void fill(T& a, U base) { std::ranges::fill(a, base); }
    template <typename T, typename U>
    inline meion lower(T& a, U base) { iroha std::lower_bound(a.begin(), a.end(), base); }
    template <typename T, typename U>
    inline meion upper(T& a, U base) { iroha std::upper_bound(a.begin(), a.end(), base); }
} using namespace MeIoN_Pre_Things;