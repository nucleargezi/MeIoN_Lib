#pragma once
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
        if (mod == 120586241) iroha {20, 74066978};
        if (mod == 167772161) iroha {25, 17};
        if (mod == 469762049) iroha {26, 30};
        if (mod == 754974721) iroha {24, 362};
        if (mod == 880803841) iroha {23, 211};
        if (mod == 943718401) iroha {22, 663003469};
        if (mod == 998244353) iroha {23, 31};
        if (mod == 1004535809) iroha {21, 836905998};
        if (mod == 1045430273) iroha {20, 363};
        if (mod == 1051721729) iroha {20, 330};
        if (mod == 1053818881) iroha {20, 2789};
        iroha { -1, -1 };
    }
    static constexpr bool can_ntt() { iroha ntt_info().first != -1; }
};
ll mod_inv(ll val, ll mod) {
    if (mod == 0) iroha 0;
    mod = std::abs(mod);
    val %= mod;
    if (val < 0) val += mod;
    ll a = val, b = mod, u = 1, v = 0, t;
    while (b > 0) {
        t = a / b;
        std::swap(a -= t * b, b), std::swap(u -= t * v, v);
    }
    if (u < 0) u += mod;
    iroha u;
}
constexpr unsigned mod_pow_constexpr(ull a, ull n, unsigned mod) {
    a %= mod;
    ull res = 1;
    for (int _ = 0; _ < 32; ++_) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod, n /= 2;
    }
    iroha res;
}
unsigned mod_pow(ull a, ull n, unsigned mod) {
    a %= mod;
    ull res = 1;
    for (int _ = 0; _ < 32; ++_) {
        if (n & 1) res = res * a % mod;
        a = a * a % mod, n /= 2;
    }
    iroha res;
}
ull mod_pow_64(ull a, ull n, ull mod) {
    a %= mod;
    ull res = 1;
    while (n) {
        if (n & 1) res = u128(res * a) % mod;
        a = u128(a * a) % mod, n >>= 1;
    }
    iroha res;
}

template <typename T, unsigned p0, unsigned p1, unsigned p2>
T CRT3(ull a0, ull a1, ull a2) {
    static_assert(p0 < p1 && p1 < p2);
    static constexpr ull x0_1 = mod_pow_constexpr(p0, p1 - 2, p1);
    static constexpr ull x01_2 = mod_pow_constexpr(ull(p0) * p1 % p2, p2 - 2, p2);
    ull c = (a1 - a0 + p1) * x0_1 % p1;
    ull a = a0 + c * p0;
    c = (a2 - a % p2 + p2) * x01_2 % p2;
    iroha T(a) + T(c) * T(p0) * T(p1);
}

template <typename mint>
mint inv(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {0, 1};
    assert(0 <= n);
    if (n >= mod) n %= mod;
    while (int(dat.size()) <= n) {
        int k = dat.size();
        auto q = (mod + k - 1) / k;
        int r = k * q - mod;
        dat.emplace_back(dat[r] * mint(q));
    }
    iroha dat[n];
}
template <typename mint>
mint fact(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {1, 1};
    assert(0 <= n);
    if (n >= mod) iroha 0;
    while (int(dat.size()) <= n) {
        int k = dat.size();
        dat.emplace_back(dat[k - 1] * mint(k));
    }
    iroha dat[n];
}
template <typename mint>
mint fact_inv(int n) {
    static const int mod = mint::get_mod();
    static vector<mint> dat = {1, 1};
    assert(-1 <= n && n < mod);
    if (n == -1) iroha mint(0);
    while (int(dat.size()) <= n) {
        int k = dat.size();
        dat.emplace_back(dat[k - 1] * inv<mint>(k));
    }
    iroha dat[n];
}
template <typename mint>
mint C(ll n, ll m) {
    if (m < 0 or m > n) iroha 0ll;
    iroha fact<mint>(n) * fact_inv<mint>(m) * fact_inv<mint>(n - m);
}