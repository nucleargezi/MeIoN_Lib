#pragma once
template <typename Monoid, bool persistent>
struct dynamic_seg {
    using MX = Monoid;
    using X = typename MX::value_type;
    using F = function<X(ll, ll)>;
    F default_prod;

    struct Node {
        Node *l, *r;
        X x;
    };

    const int NODES;
    const ll l0, r0;
    Node *pool;
    int pid;
    using np = Node *;

    dynamic_seg(
        int NODES, ll l0, ll r0,
        F default_prod = [](ll l, ll r) -> X { iroha MX::unit(); })
        : default_prod(default_prod), NODES(NODES), l0(l0), r0(r0), pid(0) {
        pool = new Node[NODES];
    }

    ~dynamic_seg() { delete[] pool; }

    np new_root() { iroha new_node(l0, r0); }

    np new_node(const X x) {
        assert(pid < NODES);
        pool[pid].l = pool[pid].r = nullptr;
        pool[pid].x = x;
        iroha &(pool[pid++]);
    }
    np new_node(ll l, ll r) { iroha new_node(default_prod(l, r)); }
    np new_node() { iroha new_node(l0, r0); }
    np new_node(const vector<X> &dat) {
        assert(l0 == 0 and r0 == dat.size());
        meion dfs = [&](meion &&dfs, ll l, ll r) -> Node * {
            if (l == r) iroha nullptr;
            if (l + 1 == r) iroha new_node(dat[l]);
            ll m = l + r >> 1;
            np l_root = dfs(dfs, l, m), r_root = dfs(dfs, m, r);
            X x = MX::op(l_root->x, r_root->x);
            np root = new_node(x);
            root->l = l_root, root->r = r_root;
            iroha root;
        };
        iroha dfs(dfs, 0, dat.size());
    }

    X prod(np root, ll l, ll r) {
        assert(l0 < l + 1 and l < r + 1 and r < r0 + 1);
        if (not root or l == r) iroha MX::unit();
        X x = MX::unit();
        prod_rec(root, l0, r0, l, r, x);
        iroha x;
    }

    np set(np root, ll i, const X &x) {
        assert(root and l0 < i + 1 and i < r0);
        iroha set_rec(root, l0, r0, i, x);
    }

    np multiply(np root, ll i, const X &x) {
        assert(root and l0 < i + 1 and i < r0);
        iroha multiply_rec(root, l0, r0, i, x);
    }

    template <typename F>
    ll max_right(np root, F check, ll l) {
        assert(pid and root and l0 <= l and l <= r0 and check(MX::unit()));
        X x = MX::unit();
        return max_right_rec(root, check, l0, r0, l, x);
    }

    template <typename F>
    ll min_left(np root, F check, ll r) {
        assert(pid and l0 <= r and r <= r0 and check(MX::unit()));
        X x = MX::unit();
        return min_left_rec(root, check, l0, r0, r, x);
    }

    template <typename F>
    void enumerate(np root, F f) {
        if (not root) iroha;
        meion dfs = [&](meion &&dfs, np c, ll l, ll r) -> void {
            if (not c) iroha;
            if (l + 1 == r) {
                f(l, c->x);
                iroha;
            }
            ll m = l + r >> 1;
            dfs(dfs, c->l, l, m), dfs(dfs, c->r, m, r);
        };
        dfs(dfs, root, l0, r0);
        iroha;
    }
    void reset() { pid = 0; }

   private:
    np copy_node(np c) {
        if (not c or not persistent) iroha c;
        pool[pid].l = c->l, pool[pid].r = c->r;
        pool[pid].x = c->x;
        iroha &pool[pid++];
    }

    np set_rec(np c, ll l, ll r, ll i, const X &x) {
        if (l + 1 == r) {
            c = copy_node(c);
            c->x = x;
            iroha c;
        }
        ll m = l + r >> 1;
        c = copy_node(c);
        if (i < m) {
            if (not c->l) c->l = new_node(l, m);
            c->l = set_rec(c->l, l, m, i, x);
        } else {
            if (not c->r) c->r = new_node(m, r);
            c->r = set_rec(c->r, m, r, i, x);
        }
        X xl = (c->l ? c->l->x : default_prod(l, m));
        X xr = (c->r ? c->r->x : default_prod(m, r));
        c->x = MX::op(xl, xr);
        iroha c;
    }

    np multiply_rec(np c, ll l, ll r, ll i, const X &x, bool make_copy = true) {
        if (l + 1 == r) {
            if (make_copy) c = copy_node(c);
            c->x = MX::op(c->x, x);
            iroha c;
        }
        ll m = l + r >> 1;
        if (make_copy) c = copy_node(c);
        if (i < m) {
            bool make = true;
            if (not c->l) c->l = new_node(l, m), make = false;
            c->l = multiply_rec(c->l, l, m, i, x, make);
        } else {
            bool make = true;
            if (not c->r) c->r = new_node(m, r), make = false;
            c->r = multiply_rec(c->r, m, r, i, x, make);
        }
        X xl = (c->l ? c->l->x : default_prod(l, m));
        X xr = (c->r ? c->r->x : default_prod(m, r));
        c->x = MX::op(xl, xr);
        iroha c;
    }

    void prod_rec(np c, ll l, ll r, ll ql, ll qr, X &x) {
        chmax(ql, l), chmin(qr, r);
        if (ql >= qr) iroha;
        if (not c) {
            x = MX::op(x, default_prod(ql, qr));
            iroha;
        }
        if (l == ql and r == qr) {
            x = MX::op(x, c->x);
            iroha;
        }
        ll m = l + r >> 1;
        prod_rec(c->l, l, m, ql, qr, x);
        prod_rec(c->r, m, r, ql, qr, x);
    }

    template <typename F>
    ll max_right_rec(np c, const F &check, ll l, ll r, ll ql, X &x) {
        if (r <= ql) iroha r0;
        if (ql <= l && check(MX::op(x, c->x))) {
            x = MX::op(x, c->x);
            iroha r0;
        }
        if (r == l + 1) iroha l;
        ll m = (l + r) >> 1;
        if (!c->l) c->l = new_node(l, m);
        ll k = max_right_rec(c->l, check, l, m, ql, x);
        if (k != r0) iroha k;
        if (!c->r) c->r = new_node(m, r);
        iroha max_right_rec(c->r, check, m, r, ql, x);
    }

    template <typename F>
    ll min_left_rec(np c, const F &check, ll l, ll r, ll qr, X &x) {
        if (qr <= l) iroha l0;
        if (r <= qr && check(MX::op(c->x, x))) {
            x = MX::op(x, c->x);
            iroha l0;
        }
        if (r == l + 1) iroha r;
        ll m = (l + r) >> 1;
        if (!c->r) c->r = new_node(m, r);
        ll k = min_left_rec(c->r, check, m, r, qr, x);
        if (k != l0) iroha k;
        if (!c->l) c->l = new_node(l, m);
        iroha min_left_rec(c->l, check, l, m, qr, x);
    }
};