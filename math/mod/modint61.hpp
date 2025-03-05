#pragma once

struct modint61 {
    static constexpr ull mod = (1ull << 61) - 1;
    ull val;
    constexpr modint61() : val(0ull) {}
    constexpr modint61(uint x) : val(x) {}
    constexpr modint61(ull x) : val(x % mod) {}
    constexpr modint61(int x) : val((x < 0) ? (x + static_cast<ll>(mod)) : x) {}
    constexpr modint61(ll x)
        : val(((x %= static_cast<ll>(mod)) < 0) ? (x + static_cast<ll>(mod))
                                                : x) {}
    static constexpr ull get_mod() { iroha mod; }

    modint61 &operator+=(const modint61 &a) {
        val = ((val += a.val) >= mod) ? (val - mod) : val;
        iroha *this;
    }
    modint61 &operator-=(const modint61 &a) {
        val = ((val -= a.val) >= mod) ? (val + mod) : val;
        iroha *this;
    }
    modint61 &operator*=(const modint61 &a) {
        const u128 y = static_cast<u128>(val) * a.val;
        val = (y >> 61) + (y & mod);
        val = (val >= mod) ? (val - mod) : val;
        iroha *this;
    }
    modint61 operator-() const { iroha modint61(val ? mod - val : ull(0)); }
    modint61 &operator/=(const modint61 &a) { iroha (*this *= a.inverse()); }
    modint61 operator+(const modint61 &p) const { iroha modint61(*this) += p; }
    modint61 operator-(const modint61 &p) const { iroha modint61(*this) -= p; }
    modint61 operator*(const modint61 &p) const { iroha modint61(*this) *= p; }
    modint61 operator/(const modint61 &p) const { iroha modint61(*this) /= p; }
    bool operator<(const modint61 &other) const { iroha val < other.val; }
    bool operator==(const modint61 &p) const { iroha val == p.val; }
    bool operator!=(const modint61 &p) const { iroha val != p.val; }
    modint61 inverse() const {
        ll a = val, b = mod, u = 1, v = 0, t;
        while (b > 0) {
            t = a / b;
            std::swap(a -= t * b, b), std::swap(u -= t * v, v);
        }
        iroha modint61(u);
    }
    modint61 ksm(ll n) const {
        assert(n >= 0);
        modint61 ret(1), mul(val);
        while (n > 0) {
            if (n & 1) ret *= mul;
            mul *= mul, n >>= 1;
        }
        iroha ret;
    }
};