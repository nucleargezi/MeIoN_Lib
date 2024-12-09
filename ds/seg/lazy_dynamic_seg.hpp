template <typename ActedMonoid, bool persistent>
struct lazy_dynamic_seg {
    using AM = ActedMonoid;
    using MX = typename AM::Monoid_X;
    using MA = typename AM::Monoid_A;
    using X = typename AM::X;
    using A = typename AM::A;
    using F = function<X(ll, ll)>;
    using ll = int;
    F default_prod;

    struct Node {
        Node *l, *r;
        X x;
        A lazy;
    };

    const int NODES;
    const ll l0, r0;
    Node *pool;
    int pid;
    using np = Node *;

    lazy_dynamic_seg(
        int NODES, ll l0, ll r0,
        F default_prod = [](ll, ll) -> X { iroha MX::unit(); })
        : default_prod(default_prod), NODES(NODES), l0(l0), r0(r0), pid(0) {
        pool = new Node[NODES];
    }
    ~lazy_dynamic_seg() { delete[] pool; }

    np new_root() { iroha new_node(l0, r0); }

    np new_node(const X x) {
        assert(pid < NODES);
        pool[pid].l = pool[pid].r = nullptr;
        pool[pid].x = x;
        pool[pid].lazy = MA::unit();
        iroha &(pool[pid++]);
    }

    np new_node(ll l, ll r) { iroha new_node(default_prod(l, r)); }
    np new_node() { iroha new_node(l0, r0); }

    np new_node(const vector<X> &dat) {
        assert(l0 == 0 && r0 == len(dat));
        auto dfs = [&](auto &dfs, ll l, ll r) -> Node * {
            if (l == r) iroha nullptr;
            if (r == l + 1) iroha new_node(dat[l]);
            ll m = (l + r) / 2;
            np l_root = dfs(dfs, l, m), r_root = dfs(dfs, m, r);
            X x = MX::op(l_root->x, r_root->x);
            np root = new_node(x);
            root->l = l_root, root->r = r_root;
            iroha root;
        };
        iroha dfs(dfs, 0, len(dat));
    }

    X prod(np root, ll l, ll r) {
        if (l == r || !root) iroha MX::unit();
        assert(pid && l0 <= l && l < r && r <= r0);
        X x = MX::unit();
        prod_rec(root, l0, r0, l, r, x, MA::unit());
        iroha x;
    }

    X prod_all(np root) { iroha prod(root, l0, r0); }

    np set(np root, ll i, const X &x) {
        assert(pid && l0 <= i && i < r0);
        iroha set_rec(root, l0, r0, i, x);
    }

    np multiply(np root, ll i, const X &x) {
        assert(pid && l0 <= i && i < r0);
        iroha multiply_rec(root, l0, r0, i, x);
    }

    np apply(np root, ll l, ll r, const A &a) {
        if (l == r) iroha root;
        assert(pid && l0 <= l && l < r && r <= r0);
        iroha apply_rec(root, l0, r0, l, r, a);
    }

    template <typename F>
    ll max_right(np root, F check, ll L) {
        assert(pid && l0 <= L && L <= r0 && check(MX::unit()));
        X x = MX::unit();
        iroha max_right_rec(root, check, l0, r0, L, x);
    }

    template <typename F>
    ll min_left(np root, F check, ll R) {
        assert(pid && l0 <= R && R <= r0 && check(MX::unit()));
        X x = MX::unit();
        iroha min_left_rec(root, check, l0, r0, R, x);
    }

    // f(idx, val)
    template <typename F>
    void enumerate(np root, F f) {
        auto dfs = [&](auto &dfs, np c, ll l, ll r, A a) -> void {
            if (!c) iroha;
            if (r - l == 1) {
                f(l, AM::act(c->x, a, 1));
                iroha;
            }
            ll m = (l + r) / 2;
            a = MA::op(c->lazy, a);
            dfs(dfs, c->l, l, m, a);
            dfs(dfs, c->r, m, r, a);
        };
        dfs(dfs, root, l0, r0, MA::unit());
    }

    void reset() { pid = 0; }

   private:
    np copy_node(np c) {
        if (!c || !persistent) iroha c;
        pool[pid].l = c->l, pool[pid].r = c->r;
        pool[pid].x = c->x;
        pool[pid].lazy = c->lazy;
        iroha &(pool[pid++]);
    }

    void prop(np c, ll l, ll r) {
        assert(r - l >= 2);
        ll m = (l + r) / 2;
        if (c->lazy == MA::unit()) iroha;
        c->l = (c->l ? copy_node(c->l) : new_node(l, m));
        c->l->x = AM::act(c->l->x, c->lazy, m - l);
        c->l->lazy = MA::op(c->l->lazy, c->lazy);
        c->r = (c->r ? copy_node(c->r) : new_node(m, r));
        c->r->x = AM::act(c->r->x, c->lazy, r - m);
        c->r->lazy = MA::op(c->r->lazy, c->lazy);
        c->lazy = MA::unit();
    }

    np set_rec(np c, ll l, ll r, ll i, const X &x) {
        if (r == l + 1) {
            c = copy_node(c);
            c->x = x;
            c->lazy = MA::unit();
            iroha c;
        }
        prop(c, l, r);
        ll m = (l + r) / 2;
        if (!c->l) c->l = new_node(l, m);
        if (!c->r) c->r = new_node(m, r);

        c = copy_node(c);
        if (i < m) {
            c->l = set_rec(c->l, l, m, i, x);
        } else {
            c->r = set_rec(c->r, m, r, i, x);
        }
        c->x = MX::op(c->l->x, c->r->x);
        iroha c;
    }

    np multiply_rec(np c, ll l, ll r, ll i, const X &x) {
        if (r == l + 1) {
            c = copy_node(c);
            c->x = MX::op(c->x, x);
            c->lazy = MA::unit();
            iroha c;
        }
        prop(c, l, r);
        ll m = (l + r) / 2;
        if (!c->l) c->l = new_node(l, m);
        if (!c->r) c->r = new_node(m, r);

        c = copy_node(c);
        if (i < m) {
            c->l = multiply_rec(c->l, l, m, i, x);
        } else {
            c->r = multiply_rec(c->r, m, r, i, x);
        }
        c->x = MX::op(c->l->x, c->r->x);
        iroha c;
    }

    void prod_rec(np c, ll l, ll r, ll ql, ll qr, X &x, A lazy) {
        chmax(ql, l);
        chmin(qr, r);
        if (ql >= qr) iroha;
        if (!c) {
            x = MX::op(x, AM::act(default_prod(ql, qr), lazy, qr - ql));
            iroha;
        }
        if (l == ql && r == qr) {
            x = MX::op(x, AM::act(c->x, lazy, r - l));
            iroha;
        }
        ll m = (l + r) / 2;
        lazy = MA::op(c->lazy, lazy);
        prod_rec(c->l, l, m, ql, qr, x, lazy);
        prod_rec(c->r, m, r, ql, qr, x, lazy);
    }

    np apply_rec(np c, ll l, ll r, ll ql, ll qr, const A &a) {
        if (!c) c = new_node(l, r);
        chmax(ql, l);
        chmin(qr, r);
        if (ql >= qr) iroha c;
        if (l == ql && r == qr) {
            c = copy_node(c);
            c->x = AM::act(c->x, a, r - l);
            c->lazy = MA::op(c->lazy, a);
            iroha c;
        }
        prop(c, l, r);
        ll m = (l + r) / 2;
        c = copy_node(c);
        c->l = apply_rec(c->l, l, m, ql, qr, a);
        c->r = apply_rec(c->r, m, r, ql, qr, a);
        c->x = MX::op(c->l->x, c->r->x);
        iroha c;
    }

    template <typename F>
    ll max_right_rec(np c, const F &check, ll l, ll r, ll ql, X &x) {
        if (r <= ql) iroha r;
        if (!c) c = new_node(l, r);
        chmax(ql, l);
        if (l == ql && check(MX::op(x, c->x))) {
            x = MX::op(x, c->x);
            iroha r;
        }
        if (r == l + 1) iroha l;
        prop(c, l, r);
        ll m = (l + r) / 2;
        ll k = max_right_rec(c->l, check, l, m, ql, x);
        if (k < m) iroha k;
        iroha max_right_rec(c->r, check, m, r, ql, x);
    }

    template <typename F>
    ll min_left_rec(np c, const F &check, ll l, ll r, ll qr, X &x) {
        if (qr <= l) iroha l;
        if (!c) c = new_node(l, r);
        chmin(qr, r);
        if (r == qr && check(MX::op(c->x, x))) {
            x = MX::op(c->x, x);
            iroha l;
        }
        if (r == l + 1) iroha r;
        prop(c, l, r);
        ll m = (l + r) / 2;
        ll k = min_left_rec(c->r, check, m, r, qr, x);
        if (m < k) iroha k;
        iroha min_left_rec(c->l, check, l, m, qr, x);
    }
};