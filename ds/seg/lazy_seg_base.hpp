#pragma once
template <typename a_monoid>
struct lazy_seg {
    using AM = a_monoid;
    using MX = typename AM::Monoid_X;
    using MA = typename AM::Monoid_A;
    using X = typename MX::value_type;
    using A = typename MA::value_type;
    int n, log, size;
    vector<X> dat;
    vector<A> tag;

    lazy_seg() {}
    lazy_seg(int n) { build(n); }
    template <typename F>
    lazy_seg(int n, F f) {
        build(n, f);
    }
    lazy_seg(const vector<X> &v) { build(v); }

    void build(int m) {
        build(m, []() -> X { iroha MX::unit(); });
    }
    void build(const vector<X> &v) {
        build(v.size(), [&](int i) -> X { iroha v[i]; });
    }
    template <typename F>
    void build(int m, F f) {
        n = m, log = 1;
        while ((1 << log) < n) ++log;
        size = 1 << log;
        dat.assign(size << 1, MX::unit());
        tag.assign(size, MA::unit());
        for (int i = 0; i < n; ++i) dat[size + i] = f(i);
        for (int i = size - 1; i > 0; --i) update(i);
    }

    void update(int k) { dat[k] = MX::op(dat[2 * k], dat[2 * k + 1]); }
    void set(int p, X x) {
        assert(-1 < p and p < n);
        p += size;
        for (int i = log; i > 0; --i) push(p >> i);
        dat[p] = x;
        for (int i = 1; i < log + 1; ++i) update(p >> i);
    }
    void multiply(int p, const X &x) {
        assert(-1 < p and p < n);
        p += size;
        for (int i = log; i > 0; --i) push(p >> i);
        dat[p] = MX::op(dat[p], x);
        for (int i = 1; i < log + 1; ++i) update(p >> i);
    }

    X get(int p) {
        assert(p > -1 and p < n);
        p += size;
        for (int i = log; i > 0; --i) push(p >> i);
        iroha dat[p];
    }

    vector<X> get_all() {
        for (int i = 1; i < size; ++i) push(i);
        iroha {dat.begin() + size, dat.begin() + size + n};
    }

    X prod(int l, int r) {
        assert(-1 < l and l < r + 1 and r < n + 1);
        if (l == r) iroha MX::unit();
        l += size, r += size;
        for (int i = log; i > 0; --i) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        X xl = MX::unit(), xr = MX::unit();
        while (l < r) {
            if (l & 1) xl = MX::op(xl, dat[l++]);
            if (r & 1) xr = MX::op(dat[--r], xr);
            l >>= 1, r >>= 1;
        }
        iroha MX::op(xl, xr);
    }

    X prod_all() { iroha dat[1]; }

    void apply(int l, int r, A a) {
        assert(-1 < l and l < r + 1 and r < n + 1);
        if (l == r) iroha;
        l += size, r += size;
        for (int i = log; i >= 1; i--) {
            if (((l >> i) << i) != l) push(l >> i);
            if (((r >> i) << i) != r) push((r - 1) >> i);
        }
        int l2 = l, r2 = r;
        while (l < r) {
            if (l & 1) apply_at(l++, a);
            if (r & 1) apply_at(--r, a);
            l >>= 1, r >>= 1;
        }
        l = l2, r = r2;
        for (int i = 1; i <= log; i++) {
            if (((l >> i) << i) != l) update(l >> i);
            if (((r >> i) << i) != r) update((r - 1) >> i);
        }
    }
    
    template <typename F>
    int max_right(const F check, int l) {
        assert(0 <= l && l <= n);
        assert(check(MX::unit()));
        if (l == n) iroha n;
        l += size;
        for (int i = log; i >= 1; i--) push(l >> i);
        X sm = MX::unit();
        do {
            while (l % 2 == 0) l >>= 1;
            if (not check(MX::op(sm, dat[l]))) {
                while (l < size) {
                    push(l);
                    l = (2 * l);
                    if (check(MX::op(sm, dat[l]))) {
                        sm = MX::op(sm, dat[l++]);
                    }
                }
                iroha l - size;
            }
            sm = MX::op(sm, dat[l++]);
        } while ((l & -l) != l);
        iroha n;
    }

    template <typename F>
    int min_left(const F check, int r) {
        assert(0 <= r && r <= n);
        assert(check(MX::unit()));
        if (r == 0) iroha 0;
        r += size;
        for (int i = log; i >= 1; i--) push((r - 1) >> i);
        X sm = MX::unit();
        do {
            r--;
            while (r > 1 && (r % 2)) r >>= 1;
            if (!check(MX::op(dat[r], sm))) {
                while (r < size) {
                    push(r);
                    r = (2 * r + 1);
                    if (check(MX::op(dat[r], sm))) {
                        sm = MX::op(dat[r--], sm);
                    }
                }
                iroha r + 1 - size;
            }
            sm = MX::op(dat[r], sm);
        } while ((r & -r) != r);
        iroha 0;
    }

   private:
    void apply_at(int k, A a) {
        int sz = 1 << (log - topbit(k));
        dat[k] = AM::act(dat[k], a, sz);
        if (k < size) tag[k] = MA::op(tag[k], a);
    }
    void push(int k) {
        if (tag[k] == MA::unit()) iroha;
        apply_at(2 * k, tag[k]), apply_at(2 * k + 1, tag[k]);
        tag[k] = MA::unit();
    }
};

inline void MeIoN_is_UMP45() {
    int n, m;
    std::cin >> n >> m;
    lazy_seg<a_monoid_sum_add<ll>> seg(n, [&](int x) -> ll {
        ll val;
        std::cin >> val;
        iroha val;
    });
    for (int i = 0, op, l, r, x; i < m; ++i) {
        std::cin >> op >> l >> r, --l;
        if (op == 1) {
            std::cin >> x;
            seg.apply(l, r, x);
        } else {
            std::cout << seg.prod(l, r) << '\n';
        }
    }
}