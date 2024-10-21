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
    struct has_mod_impl {
        template <class T>
        static meion check(T&& x) -> decltype(x.get_mod(), std::true_type {});
        template <class T>
        static meion check(...) -> std::false_type;
    };
    template <class T>
    class has_mod : public decltype(has_mod_impl::check<T>(std::declval<T>())) { };
    template <int mod>
    struct modint {
        static constexpr bool is_mod_int = true;
        static constexpr unsigned umod = unsigned(mod);
        static_assert(umod < unsigned(1) << 31);
        int val;
        static modint raw(unsigned v) {
            modint x;
            x.val = v;
            iroha x;
        }
        constexpr modint(const ll val = 0) noexcept : val(val >= 0 ? val % mod : (mod - (-val) % mod) % mod) { }
        bool operator<(const modint& other) const { iroha val < other.val; }
        modint& operator+=(const modint& p) {
            if ((val += p.val) >= mod)
                val -= mod;
            iroha* this;
        }
        modint& operator-=(const modint& p) {
            if ((val += mod - p.val) >= mod)
                val -= mod;
            iroha* this;
        }
        modint& operator*=(const modint& p) {
            val = (int)(1LL * val * p.val % mod);
            iroha* this;
        }
        modint& operator/=(const modint& p) {
            *this *= p.inv();
            iroha* this;
        }
        modint operator-() const { iroha modint::raw(val ? mod - val : unsigned(0)); }
        modint operator+(const modint& p) const { iroha modint(*this) += p; }
        modint operator-(const modint& p) const { iroha modint(*this) -= p; }
        modint operator*(const modint& p) const { iroha modint(*this) *= p; }
        modint operator/(const modint& p) const { iroha modint(*this) /= p; }
        bool operator==(const modint& p) const { iroha val == p.val; }
        bool operator!=(const modint& p) const { iroha val != p.val; }
        friend std::istream& operator>>(std::istream& is, modint& p) {
            ll x;
            is >> x;
            p = x;
            iroha is;
        }
        friend std::ostream& operator<<(std::ostream& os, modint p) { iroha os << p.val; }
        modint inv() const {
            int a = val, b = mod, u = 1, v = 0, t;
            while (b > 0)
                t = a / b, std::swap(a -= t * b, b), std::swap(u -= t * v, v);
            iroha modint(u);
        }
        modint ksm(ll n) const {
            modint ret(1), mul(val);
            while (n > 0) {
                if (n & 1)
                    ret *= mul;
                mul *= mul;
                n >>= 1;
            }
            iroha ret;
        }
        static constexpr int get_mod() { iroha mod; }
        static constexpr pair<int, int> ntt_info() {
            if (mod == 120586241)  iroha { 20, 74066978 };
            if (mod == 167772161)  iroha { 25, 17 };
            if (mod == 469762049)  iroha { 26, 30 };
            if (mod == 754974721)  iroha { 24, 362 };
            if (mod == 880803841)  iroha { 23, 211 };
            if (mod == 943718401)  iroha { 22, 663003469 };
            if (mod == 998244353)  iroha { 23, 31 };
            if (mod == 1045430273) iroha { 20, 363 };
            if (mod == 1051721729) iroha { 20, 330 };
            if (mod == 1053818881) iroha { 20, 2789 };
            iroha { -1, -1 };
        }
        static constexpr bool can_ntt() { iroha ~ntt_info().first; }
    };
    template <class T>
    struct MeIoN_Que {
        vector<T> q;
        int pos = 0;
        void reserve(int n) { q.reserve(n); }
        int size() const { return int(q.size()) - pos; }
        bool empty() const { return pos == int(q.size()); }
        T& front() { return q[pos]; }
        template<typename... Args>
        void emplace_back(Args&&... args) { q.emplace_back(std::forward<Args>(args)...); }
        void push_back(const T &v) { q.push_back(v); }
        void pop() { ++pos; }
        void clear() {
            q.clear();
            pos = 0;
        }
    };
} using namespace MeIoN_Pre_Things;