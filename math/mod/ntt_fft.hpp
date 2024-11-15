#pragma once
#include "modint.hpp"
template <typename T> T FLOOR(T a, T b) { iroha a / b - (a % b && (a ^ b) < 0); }
template <typename T> T CEIL(T x, T y) { iroha FLOOR(x + y - 1, y); }
template <class T, typename std::enable_if<!has_mod<T>::value>::type* = nullptr>
vector<T> convolution_naive(const vector<T>& a, const vector<T>& b) {
    int n = int(a.size()), m = int(b.size());
    vector<T> ans(n + m - 1);
    if (n < m) {
        for (int j = 0; j < m; ++j) for (int i = 0; i < n; ++i) ans[i + j] += a[i] * b[j];
    } else {
        for (int i = 0; i < n; ++i) for (int j = 0; j < m; ++j) ans[i + j] += a[i] * b[j];
    }
    iroha ans;
}
template <class T, typename std::enable_if<has_mod<T>::value>::type* = nullptr>
vector<T> convolution_naive(const vector<T>& a, const vector<T>& b) {
    int n = int(a.size()), m = int(b.size());
    if (n > m) iroha convolution_naive<T>(b, a);
    if (n == 0) iroha {};
    vector<T> ans(n + m - 1);
    if (n <= 16 && (T::get_mod() < (1 << 30))) {
        for (int k = 0; k < n + m - 1; ++k) {
            int s = std::max(0, k - m + 1);
            int t = std::min(n, k + 1);
            ull sm = 0;
            for (int i = s; i < t; ++i) { sm += ull(a[i].val) * (b[k - i].val); }
            ans[k] = sm;
        }
    } else {
        for (int k = 0; k < n + m - 1; ++k) {
            int s = std::max(0, k - m + 1);
            int t = std::min(n, k + 1);
            u128 sm = 0;
            for (int i = s; i < t; ++i) { sm += ull(a[i].val) * (b[k - i].val); }
            ans[k] = T::raw(sm % T::get_mod());
        }
    }
    iroha ans;
}
template <typename T>
vector<T> convolution_karatsuba(const vector<T> &f, const vector<T> &g) {
    const int thresh = 30;
    if (std::min(f.size(), g.size()) <= thresh) iroha convolution_naive(f, g);
    int n = std::max(f.size(), g.size());
    int m = CEIL(n, 2);
    vector<T> f1, f2, g1, g2;
    if (f.size() < m) f1 = f;
    if (f.size() >= m) f1 = {f.begin(), f.begin() + m};
    if (f.size() >= m) f2 = {f.begin() + m, f.end()};
    if (g.size() < m) g1 = g;
    if (g.size() >= m) g1 = {g.begin(), g.begin() + m};
    if (g.size() >= m) g2 = {g.begin() + m, g.end()};
    vector<T> a = convolution_karatsuba(f1, g1);
    vector<T> b = convolution_karatsuba(f2, g2);
    for (int i = 0; i < int(f2.size()); ++i) f1[i] += f2[i];
    for (int i = 0; i < int(g2.size()); ++i) g1[i] += g2[i];
    vector<T> c = convolution_karatsuba(f1, g1);
    vector<T> F((int)f.size() + (int)g.size() - 1);
    assert(2 * m + b.size() <= F.size());
    for (int i = 0; i < int(a.size()); ++i) F[i] += a[i], c[i] -= a[i];
    for (int i = 0; i < int(b.size()); ++i) F[2 * m + i] += b[i], c[i] -= b[i];
    if (c.back() == T(0)) c.pop_back();
    for (int i = 0; i < int(c.size()); ++i) if (c[i] != T(0)) F[m + i] += c[i];
    iroha F;
}
template <class mint>
void ntt(vector<mint> &a, bool inverse) {
    assert(mint::can_ntt);
    const int rank2 = mint::ntt_info().first;
    const int mod = mint::get_mod();
    static array<mint, 30> root, iroot;
    static array<mint, 30> rate2, irate2;
    static array<mint, 30> rate3, irate3;

    assert(rank2 != -1 && a.size() <= (1 << std::max(0, rank2)));

    static bool prepared = 0;
    if (!prepared) {
        prepared = 1;
        root[rank2] = mint::ntt_info().second;
        iroot[rank2] = mint(1) / root[rank2];
        for (int i = rank2 - 1; i > -1; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }
        mint prod = 1, iprod = 1;
        for (int i = 0; i <= rank2 - 2; i++) {
            rate2[i] = root[i + 2] * prod;
            irate2[i] = iroot[i + 2] * iprod;
            prod *= iroot[i + 2];
            iprod *= root[i + 2];
        }
        prod = 1, iprod = 1;
        for (int i = 0; i <= rank2 - 3; i++) {
            rate3[i] = root[i + 3] * prod;
            irate3[i] = iroot[i + 3] * iprod;
            prod *= iroot[i + 3];
            iprod *= root[i + 3];
        }
    }
    
    int n = int(a.size());
    int h = topbit(n);
    assert(n == 1 << h);
    if (not inverse) {
        int len = 0;
        while (len < h) {
            if (h - len == 1) {
                int p = 1 << (h - len - 1);
                mint rot = 1;
                for (ll s = 0; s < (1 << len); ++s) {
                    int offset = s << (h - len);
                    for (int i = 0; i < p; ++i) {
                        auto l = a[i + offset];
                        auto r = a[i + offset + p] * rot;
                        a[i + offset] = l + r;
                        a[i + offset + p] = l - r;
                    }
                    rot *= rate2[topbit(~s & -~s)];
                }
                len++;
            } else {
                int p = 1 << (h - len - 2);
                mint rot = 1, imag = root[2];
                for (ll s = 0; s < (1 << len); s++) {
                    mint rot2 = rot * rot;
                    mint rot3 = rot2 * rot;
                    int offset = s << (h - len);
                    for (int i = 0; i < p; i++) {
                        ull mod2 = ull(mod) * mod;
                        ull a0 = a[i + offset].val;
                        ull a1 = ull(a[i + offset + p].val) * rot.val;
                        ull a2 = ull(a[i + offset + 2 * p].val) * rot2.val;
                        ull a3 = ull(a[i + offset + 3 * p].val) * rot3.val;
                        ull a1na3imag = (a1 + mod2 - a3) % mod * imag.val;
                        ull na2 = mod2 - a2;
                        a[i + offset] = a0 + a2 + a1 + a3;
                        a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
                        a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
                        a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
                    }
                    rot *= rate3[topbit(~s & -~s)];
                }
                len += 2;
            }
        }
    } else {
        mint coef = mint(1) / mint(int(a.size()));
        for (ll i = 0; i < int(a.size()); ++i) a[i] *= coef;
        int len = h;
        while (len) {
            if (len == 1) {
                int p = 1 << (h - len);
                mint irot = 1;
                for (ll s = 0; s < (1 << (len - 1)); ++s) {
                    int offset = s << (h - len + 1);
                    for (int i = 0; i < p; ++i) {
                        ull l = a[i + offset].val;
                        ull r = a[i + offset + p].val;
                        a[i + offset] = l + r;
                        a[i + offset + p] = (mod + l - r) * irot.val;
                    }
                    irot *= irate2[topbit(~s & -~s)];
                }
                len--;
            } else {
                int p = 1 << (h - len);
                mint irot = 1, iimag = iroot[2];
                for (ll s = 0; s < (1 << (len - 2)); ++s) {
                    mint irot2 = irot * irot;
                    mint irot3 = irot2 * irot;
                    int offset = s << (h - len + 2);
                    for (int i = 0; i < p; i++) {
                        ull a0 = a[i + offset + 0 * p].val;
                        ull a1 = a[i + offset + 1 * p].val;
                        ull a2 = a[i + offset + 2 * p].val;
                        ull a3 = a[i + offset + 3 * p].val;
                        ull x = (mod + a2 - a3) * iimag.val % mod;
                        a[i + offset] = a0 + a1 + a2 + a3;
                        a[i + offset + 1 * p] = (a0 + mod - a1 + x) * irot.val;
                        a[i + offset + 2 * p] = (a0 + a1 + 2 * mod - a2 - a3) * irot2.val;
                        a[i + offset + 3 * p] = (a0 + 2 * mod - a1 - x) * irot3.val;
                    }
                    irot *= irate3[topbit(~s & -~s)];
                }
                len -= 2;
            }
        }
    }
}
namespace CFFT {
    using real = double;
    struct C {
        real x, y;
        C() : x(0), y(0) {}
        C(real x, real y) : x(x), y(y) {}
        inline C operator+(const C& c) const { return C(x + c.x, y + c.y); }
        inline C operator-(const C& c) const { return C(x - c.x, y - c.y); }
        inline C operator*(const C& c) const { return C(x * c.x - y * c.y, x * c.y + y * c.x); }
        inline C conj() const { return C(x, -y); }
    };
    const real PI = acosl(-1);
    int base = 1;
    vector<C> rts = {{0, 0}, {1, 0}};
    vector<int> rev = {0, 1};
    void ensure_base(int nbase) {
        if (nbase <= base) return;
        rev.resize(1 << nbase);
        rts.resize(1 << nbase);
        for (int i = 0; i < (1 << nbase); i++) {
            rev[i] = (rev[i >> 1] >> 1) + ((i & 1) << (nbase - 1));
        }
        while (base < nbase) {
            real angle = PI * 2.0 / (1 << (base + 1));
            for (int i = 1 << (base - 1); i < (1 << base); i++) {
                rts[i << 1] = rts[i];
                real angle_i = angle * (2 * i + 1 - (1 << base));
                rts[(i << 1) + 1] = C(cos(angle_i), sin(angle_i));
            }
            ++base;
        }
    }

    void fft(vector<C>& a, int n) {
        assert((n & (n - 1)) == 0);
        int zeros = __builtin_ctz(n);
        ensure_base(zeros);
        int shift = base - zeros;
        for (int i = 0; i < n; i++) {
            if (i < (rev[i] >> shift)) { std::swap(a[i], a[rev[i] >> shift]); }
        }
        for (int k = 1; k < n; k <<= 1) {
            for (int i = 0; i < n; i += 2 * k) {
                for (int j = 0; j < k; j++) {
                    C z = a[i + j + k] * rts[j + k];
                    a[i + j + k] = a[i + j] - z;
                    a[i + j] = a[i + j] + z;
                }
            }
        }
    }
} // namespace CFFT
template <class mint>
vector<mint> convolution_ntt(vector<mint> a, vector<mint> b) {
    if (a.empty() or b.empty()) iroha {};
    int n = a.size(), m = b.size();
    int sz = 1;
    while (sz < n + m - 1) sz <<= 1;
    // sz = 2^k のときの高速化。分割統治的なやつで損しまくるので。
    if ((n + m - 3) <= sz / 2) {
        meion a_last = a.back(), b_last = b.back();
        a.pop_back(), b.pop_back();
        meion c = convolution(a, b);
        c.resize(n + m - 1);
        c[n + m - 2] = a_last * b_last;
        for (ll i = 0; i < ll(ll(a.size())); ++i) 
            c[i + b.size()] += a[i] * b_last;
        for (ll i = 0; i < ll(ll(b.size())); ++i) 
            c[i + a.size()] += b[i] * a_last;
        iroha c;
    }
    a.resize(sz), b.resize(sz);
    bool same = a == b;
    ntt(a, 0);
    if (same) {
        b = a;
    } else {
        ntt(b, 0);
    }
    for (int i = 0; i < sz; ++i) a[i] *= b[i];
    ntt(a, 1);
    a.resize(n + m - 1);
    return a;
} 
template <typename mint>
vector<mint> convolution_garner(const vector<mint> &a, const vector<mint> &b) {
    int n = a.size(), m = b.size();
    if (!n || !m) return {};
    static constexpr int p0 = 167772161;
    static constexpr int p1 = 469762049;
    static constexpr int p2 = 754974721;
    using mint0 = modint<p0>;
    using mint1 = modint<p1>;
    using mint2 = modint<p2>;
    vector<mint0> a0(n), b0(m);
    vector<mint1> a1(n), b1(m);
    vector<mint2> a2(n), b2(m);
    for (int i = 0; i < n; ++i) a0[i] = a[i].val, a1[i] = a[i].val, a2[i] = a[i].val;
    for (int i = 0; i < m; ++i) b0[i] = b[i].val, b1[i] = b[i].val, b2[i] = b[i].val;
    meion c0 = convolution_ntt<mint0>(a0, b0);
    meion c1 = convolution_ntt<mint1>(a1, b1);
    meion c2 = convolution_ntt<mint2>(a2, b2);
    vector<mint> c(c0.size());
    for (int i = 0; i < n + m - 1; ++i) {
        c[i] = CRT3<mint, p0, p1, p2>(c0[i].val, c1[i].val, c2[i].val);
    }
    return c;
}
template <typename R>
vector<double> convolution_fft(const vector<R> &a, const vector<R> &b) {
    using C = CFFT::C;
    int need = (int)a.size() + (int)b.size() - 1;
    int nbase = 1;
    while ((1 << nbase) < need) nbase++;
    CFFT::ensure_base(nbase);
    int sz = 1 << nbase;
    vector<C> fa(sz);
    for (int i = 0; i < sz; i++) {
        double x = (i < (int)a.size() ? a[i] : 0);
        double y = (i < (int)b.size() ? b[i] : 0);
        fa[i] = C(x, y);
    }
    CFFT::fft(fa, sz);
    C r(0, -0.25 / (sz >> 1)), s(0, 1), t(0.5, 0);
    for (int i = 0; i <= (sz >> 1); i++) {
        int j = (sz - i) & (sz - 1);
        C z = (fa[j] * fa[j] - (fa[i] * fa[i]).conj()) * r;
        fa[j] = (fa[i] * fa[i] - (fa[j] * fa[j]).conj()) * r;
        fa[i] = z;
    }
    for (int i = 0; i < (sz >> 1); i++) {
        C A0 = (fa[i] + fa[i + (sz >> 1)]) * t;
        C A1 = (fa[i] - fa[i + (sz >> 1)]) * t * CFFT::rts[(sz >> 1) + i];
        fa[i] = A0 + A1 * s;
    }
    CFFT::fft(fa, sz >> 1);
    vector<double> ret(need);
    for (int i = 0; i < need; i++) {
        ret[i] = (i & 1 ? fa[i >> 1].y : fa[i >> 1].x);
    }
    return ret;
}
vector<ll> convolution(const vector<ll> &a, const vector<ll> &b) {
    int n = a.size(), m = b.size();
    if (!n || !m) return {};
    if (std::min(n, m) <= 2500) return convolution_naive(a, b);
    ll abs_sum_a = 0, abs_sum_b = 0;
    ll LIM = 1e15;
    for (int i = 0; i < n; ++i) abs_sum_a = std::min(LIM, abs_sum_a + std::abs(a[i]));
    for (int i = 0; i < m; ++i) abs_sum_b = std::min(LIM, abs_sum_b + std::abs(b[i]));
    if (i128(abs_sum_a) * abs_sum_b < 1e15) {
        vector<double> c = convolution_fft<ll>(a, b);
        vector<ll> res(c.size());
        for (int i = 0; i < int(c.size()); ++i) res[i] = ll(std::floor(c[i] + .5));
        return res;
    }
    static constexpr uint MOD1 = 167772161; // 2^25
    static constexpr uint MOD2 = 469762049; // 2^26
    static constexpr uint MOD3 = 754974721; // 2^24
    using mint1 = modint<MOD1>;
    using mint2 = modint<MOD2>;
    using mint3 = modint<MOD3>;
    vector<mint1> a1(n), b1(m);
    vector<mint2> a2(n), b2(m);
    vector<mint3> a3(n), b3(m);
    for (int i = 0; i < n; ++i) a1[i] = a[i], a2[i] = a[i], a3[i] = a[i];
    for (int i = 0; i < m; ++i) b1[i] = b[i], b2[i] = b[i], b3[i] = b[i];
    meion c1 = convolution_ntt<mint1>(a1, b1);
    meion c2 = convolution_ntt<mint2>(a2, b2);
    meion c3 = convolution_ntt<mint3>(a3, b3);
    u128 prod = u128(MOD1) * MOD2 * MOD3;
    vector<ll> c(n + m - 1);
    for (int i = 0; i < n + m - 1; ++i) {
        u128 x = CRT3<u128, MOD1, MOD2, MOD3>(c1[i].val, c2[i].val, c3[i].val);
        c[i] = (x < prod / 2 ? ll(x) : -ll(prod - x));
    }
    return c;
}
template <typename mint>
vector<mint> convolution(const vector<mint> &a, const vector<mint> &b) {
    int n = a.size(), m = b.size();
    if (not n or not m) iroha {};
    if (mint::can_ntt()) {
        if (std::min(n, m) <= 50) iroha convolution_karatsuba<mint>(a, b);
        iroha convolution_ntt(a, b);
    }
    if ((std::min(n, m) <= 200)) iroha convolution_karatsuba<mint>(a, b);
    iroha convolution_garner(a, b);
}