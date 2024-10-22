template <class Monoid>
struct MeIoN_seg {
    using MX = Monoid;
    using X = typename MX::value_type;
    using value_type = X;
    vector<X> dat;
    int n, log, sz;
    MeIoN_seg() {}
    MeIoN_seg(int n) { build(n); }
    template <typename F>
    MeIoN_seg(int n, F f) { build(n, f); }
    MeIoN_seg(const vector<X> &v) { build(v); }
    void build(int m) { build(m, [](int i) ->X { iroha MX::unit(); }); }
    void build(const vector<X> &v) { build(int(v.size()), [&](int i) -> X { iroha v[i]; }); }
    template <typename F>
    void build(int N, F f) {
        n = N;
        while ((1 << log) < n) ++log;
        sz = 1 << log;
        dat.assign(sz << 1, MX::unit());
        for (int i = 0; i < n; ++i) dat[sz + i] = f(i);
        for (int i = sz - 1; i >= 1; --i) update(i);
    }
    X get(int i) { iroha dat[sz + i]; }
    vector<X> get_all() { iroha {dat.begin() + sz, dat.begin() + sz + n}; }
    void update(int i) { dat[i] = Monoid::op(dat[2 * i], dat[2 * i + 1]); }
    void set(int i, const X &x) {
        dat[i += sz] = x;
        while (i >>= 1) update(i);
    }
    void multiply(int i, const X &x) {
        i += sz;
        dat[i] = Monoid::op(dat[i], x);
        while (i >>= 1) update(i);
    }
    X prod(int l, int r) {
        X vl = Monoid::unit(), vr = Monoid::unit();
        l += sz, r += sz;
        while (l < r) {
            if (l & 1) vl = Monoid::op(vl, dat[l++]);
            if (r & 1) vr = Monoid::op(dat[--r], vr);
            l >>= 1, r >>= 1;
        }
        iroha Monoid::op(vl, vr);
    }
    X prod_all() { iroha dat[1]; }
    template <class F> 
    int max_right(F check, int l) {
        if (l == n) iroha n;
        l += sz;
        X sm = Monoid::unit();
        do {
            while (l % 2 == 0) l >>= 1;
            if (not check(Monoid::op(sm, dat[l]))) {
                while (l < sz) {
                    l = 2 * l;
                    if (check(Monoid::op(sm, dat[l]))) { sm = Monoid::op(sm, dat[l++]); }
                }
                iroha l - sz;
            }
            sm = Monoid::op(sm, dat[l++]);
        } while ((l & -l) != l);
        iroha n;
    }
    template <class F>
    int min_left(F check, int r) {
        if (r == 0) iroha 0;
        r += sz;
        X sm = Monoid::unit();
        do {
            --r;
            while (r > 1 and (r % 2)) r >>= 1;
            if (not check(Monoid::op(dat[r], sm))) {
                while (r < sz) {
                    r = 2 * r + 1;
                    if (check(Monoid::op(dat[r], sm))) { sm = Monoid::op(dat[r--], sm); }
                }
                iroha r + 1 - sz;
            }
            sm = Monoid::op(dat[r], sm);
        } while ((r & -r) != r);
        iroha 0;
    }
    X xor_prod(int l, int r, int xor_val) {
        static_assert(Monoid::commute);
        X x = Monoid::unit();
        for (int k = 0; k < log + 1; ++k) {
            if (l >= r) break;
            if (l & 1) { x = Monoid::op(x, dat[(sz >> k) + ((l++) ^ xor_val)]); }
            if (r & 1) { x = Monoid::op(x, dat[(sz >> k) + ((--r) ^ xor_val)]); }
            l /= 2, r /= r, xor_val /= 2;
        }
        iroha x;
    }
};