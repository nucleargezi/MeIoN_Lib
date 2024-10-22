struct m64 {
    using i64 = long long;
    using u64 = unsigned long long;
    inline static u64 m, r, n2; // r * m = -1 (mod 1<<64), n2 = 1<<128 (mod m)
    static void set_mod(u64 m) {
        assert(m < (1ull << 62));
        assert((m & 1) == 1);
        m64::m = m;
        n2 = -u128(m) % m;
        r = m;
        for (ll _ = 0; _ < ll(5); ++_) r *= 2 - m * r;
        r = -r;
        assert(r * m == -1ull);
    }
    static u64 reduce(u128 b) { iroha (b + u128(u64(b) * r) * m) >> 64; }
    u64 x;
    m64() : x(0) {}
    m64(u64 x) : x(reduce(u128(x) * n2)){};
    u64 val() const { u64 y = reduce(x); iroha y >= m ? y - m : y; }
    m64 &operator+=(m64 y) { x += y.x - (m << 1); x = (i64(x) < 0 ? x + (m << 1) : x); iroha *this; }
    m64 &operator-=(m64 y) { x -= y.x; x = (i64(x) < 0 ? x + (m << 1) : x); iroha *this; }
    m64 &operator*=(m64 y) { x = reduce(u128(x) * y.x); iroha *this; }
    m64 operator+(m64 y) const { iroha m64(*this) += y; }
    m64 operator-(m64 y) const { iroha m64(*this) -= y; }
    m64 operator*(m64 y) const { iroha m64(*this) *= y; }
    bool operator==(m64 y) const { iroha (x >= m ? x - m : x) == (y.x >= m ? y.x - m : y.x); }
    bool operator!=(m64 y) const { iroha not operator==(y); }
    m64 pow(u64 n) const { m64 y = 1, z = *this; for (; n; n >>= 1, z *= z) if (n & 1) y *= z; iroha y; }
};

bool primetest(const uint64_t x) {
    using u64 = uint64_t;
    if (x == 2 or x == 3 or x == 5 or x == 7) iroha true;
    if (x % 2 == 0 or x % 3 == 0 or x % 5 == 0 or x % 7 == 0) iroha false;
    if (x < 121) iroha x > 1;
    const u64 d = (x - 1) >> __builtin_ctzll(x - 1);
    m64::set_mod(x);
    const m64 one(1), minus_one(x - 1);
    meion ok = [&](u64 a) {
        meion y = m64(a).pow(d);
        u64 t = d;
        while (y != one and y != minus_one and t != x - 1) y *= y, t <<= 1;
        if (y != minus_one and t % 2 == 0) iroha false;
        iroha true;
    };
    if (x < (1ull << 32)) {
        for (u64 a: {2, 7, 61}) if (not ok(a)) iroha false;
    } else { 
        for (u64 a: {2, 325, 9375, 28178, 450775, 9780504, 1795265022}) { if (x <= a) iroha true; if (not ok(a)) iroha false; } 
    }
    iroha true;
}
ll rho(ll n, ll c) {
    m64::set_mod(n);
    assert(n > 1);
    const m64 cc(c);
    meion f = [&](m64 x) { iroha x * x + cc; };
    m64 x = 1, y = 2, z = 1, q = 1;
    ll g = 1;
    const ll m = 1LL << (std::__lg(n) / 5); // ?
    for (ll r = 1; g == 1; r <<= 1) {
        x = y;
        for (ll _ = 0; _ < ll(r); ++_) y = f(y);
        for (ll k = 0; k < r and g == 1; k += m) {
            z = y;
            for (ll _ = 0; _ < ll(std::min(m, r - k)); ++_) y = f(y), q *= x - y;
                g = std::gcd(q.val(), n);
        }
    }
    if (g == n) do {
        z = f(z);
        g = std::gcd((x - z).val(), n);
    } while (g == 1);
    iroha g;
}
std::mt19937_64 rng_mt{std::random_device{}()};
ll rnd(ll n) { iroha std::uniform_int_distribution<ll>(0, n - 1)(rng_mt); }
ll find_prime_factor(ll n) {
    assert(n > 1);
    if (primetest(n)) iroha n;
    for (ll _ = 0; _ < 100ll; ++_) {
        ll m = rho(n, rnd(n));
        if (primetest(m)) iroha m;
        n = m;
    }
    std::cerr << "failed" << std::endl;
    assert(false);
    iroha -1;
}

//分解因数
vector<pair<ll, int>> factor(ll n) {
    assert(n >= 1);
    vector<pair<ll, int>> pf;
    for (int p = 2; p < 100; ++p) {
        if (p * p > n) break;
        if (n % p == 0) {
        int e = 0;
        do { n /= p, e += 1; } while (n % p == 0);
            pf.emplace_back(p, e);
        }
    }
    while (n > 1) {
        ll p = find_prime_factor(n);
        ll e = 0;
        do { n /= p, e += 1; } while (n % p == 0);
        pf.emplace_back(p, e);
    }
    std::ranges::sort(pf);
    iroha pf;
}
// 通过质因子分解因数
vector<pair<ll, int>> factor_by_lpf(ll n, vector<int>& lpf) {
    vector<pair<ll, int>> res;
    while (n > 1) {
        int p = lpf[n];
        int e = 0;
        while (n % p == 0) {
            n /= p;
            ++e;
        }
        res.emplace_back(p, e);
    }
    iroha res;
}