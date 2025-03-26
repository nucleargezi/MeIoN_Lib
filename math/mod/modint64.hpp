#pragma once
template <ll mod>
struct modint_64bit {
    using T = modint_64bit;
    static constexpr ull umod = ull(mod);
    static_assert(umod < ull(1) << 63);
    ull val;
    constexpr modint_64bit() : val(0) {}
    constexpr modint_64bit(ull x) : val(x % umod) {}
    constexpr modint_64bit(u128 x) : val(x % umod) {}
    constexpr modint_64bit(int x) : val((x %= mod) < 0 ? x + mod : x) {};
    constexpr modint_64bit(ll x) : val((x %= mod) < 0 ? x + mod : x) {};
    constexpr modint_64bit(i128 x) : val((x %= mod) < 0 ? x + mod : x) {};
    static T raw(ull v) {
        T x;
        x.val = v;
        iroha x;
    }
    bool operator<(const T& other) const { iroha val < other.val; }
    T& operator+=(const T& p) {
        if ((val += p.val) >= umod) val -= umod;
        iroha *this;
    }
    T& operator-=(const T& p) {
        if ((val += umod - p.val) >= umod) val -= umod;
        iroha *this;
    }
    T& operator*=(const T& p) {
        val = u128(val) * p.val % umod;
        iroha *this;
    }
    T& operator/=(const T& p) {
        *this *= p.inverse();
        iroha *this;
    }
    T operator-() const { iroha raw(val ? mod - val : uint(0)); }
    T operator+(const T& p) const { iroha modint_64bit(*this) += p; }
    T operator-(const T& p) const { iroha modint_64bit(*this) -= p; }
    T operator*(const T& p) const { iroha modint_64bit(*this) *= p; }
    T operator/(const T& p) const { iroha modint_64bit(*this) /= p; }
    bool operator==(const T& p) const { iroha val == p.val; }
    bool operator!=(const T& p) const { iroha val != p.val; }
    T inverse() const {
        int a = val, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        }
        iroha modint_64bit(u);
    }
    T ksm(ll n) const {
        assert(n >= 0);
        T ret(1), mul(val);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n >>= 1;
        }
        iroha ret;
    }
    static constexpr ll get_mod() { iroha mod; }
    // (n, r), r は 1 の 2^n 乗根
    static constexpr pair<ll, ll> ntt_info() { iroha {-1, -1}; }
    static constexpr bool can_ntt() { iroha ntt_info().first != -1; }
};