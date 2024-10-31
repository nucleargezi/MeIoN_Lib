#pragma once

template <typename tree, typename Node, typename T, T (*f)(T, T), T (*ts)(T)>
struct rev_BBST : tree {
    using tree::merge;
    using tree::split;
    using typename tree::ptr;

    rev_BBST() = default;

    void toggle(ptr t) {
        std::swap(t->l, t->r);
        t->sum = ts(t->sum);
        t->rev ^= true;
    }

    T fold(ptr &t, int a, int b) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        auto ret = sum(y.first);
        t = merge(x.first, merge(y.first, y.second));
        return ret;
    }

    void reverse(ptr &t, int a, int b) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        toggle(y.first);
        t = merge(x.first, merge(y.first, y.second));
    }

    ptr update(ptr t) override {
        if (!t) return t;
        t->cnt = 1;
        t->sum = t->key;
        if (t->l) t->cnt += t->l->cnt, t->sum = f(t->l->sum, t->sum);
        if (t->r) t->cnt += t->r->cnt, t->sum = f(t->sum, t->r->sum);
        return t;
    }

   protected:
    inline T sum(const ptr t) { return t ? t->sum : T(); }

    void push(ptr t) override {
        if (!t) return;
        if (t->rev) {
            if (t->l) toggle(t->l);
            if (t->r) toggle(t->r);
            t->rev = false;
        }
    }
};

template <typename Tree, typename Node, typename T, typename E, T (*f)(T, T),
          T (*g)(T, E), E (*h)(E, E), T (*ts)(T)>
struct lazy_rev_BBST : Tree {
    using Tree::merge;
    using Tree::split;
    using typename Tree::ptr;

    lazy_rev_BBST() = default;

    void toggle(ptr t) {
        if (!t) return;
        swap(t->l, t->r);
        t->sum = ts(t->sum);
        t->rev ^= true;
    }

    T fold(ptr &t, int a, int b) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        auto ret = sum(y.first);
        t = merge(x.first, merge(y.first, y.second));
        return ret;
    }

    void reverse(ptr &t, int a, int b) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        toggle(y.first);
        t = merge(x.first, merge(y.first, y.second));
    }

    void apply(ptr &t, int a, int b, const E &e) {
        auto x = split(t, a);
        auto y = split(x.second, b - a);
        propagate(y.first, e);
        t = merge(x.first, merge(y.first, y.second));
    }

   protected:
    inline T sum(const ptr t) { return t ? t->sum : T(); }

    ptr update(ptr t) override {
        if (!t) return t;
        t->cnt = 1;
        t->sum = t->key;
        if (t->l) t->cnt += t->l->cnt, t->sum = f(t->l->sum, t->sum);
        if (t->r) t->cnt += t->r->cnt, t->sum = f(t->sum, t->r->sum);
        return t;
    }

    void push(ptr t) override {
        if (!t) return;
        if (t->rev) {
            if (t->l) toggle(t->l);
            if (t->r) toggle(t->r);
            t->rev = false;
        }
        if (t->lazy != E()) {
            if (t->l) propagate(t->l, t->lazy);
            if (t->r) propagate(t->r, t->lazy);
            t->lazy = E();
        }
    }

    void propagate(ptr t, const E &x) {
        t->lazy = h(t->lazy, x);
        t->key = g(t->key, x);
        t->sum = g(t->sum, x);
    }
};

template <typename T>
struct rev_splay_node {
    using ptr = rev_splay_node *;
    ptr l, r, p;
    T key, sum;
    int cnt;
    bool rev;

    rev_splay_node(const T &t = T())
        : l(), r(), p(), key(t), sum(t), cnt(1), rev(false) {}
};

template <typename T, typename E>
struct lazy_rev_splay_node {
    using ptr = lazy_rev_splay_node *;
    ptr l, r, p;
    T key, sum;
    E lazy;
    int cnt;
    bool rev;

    lazy_rev_splay_node(const T &t = T(), const E &e = E())
        : l(), r(), p(), key(t), sum(t), lazy(e), cnt(1), rev(false) {}
};

template <typename Node>
struct splay_tree_base {
    // using ptr = Node<int> *;
    using ptr = Node *;
    template <typename... Args>
    ptr my_new(const Args &...args) {
        iroha new Node(args...);
    }
    void my_del(ptr p) { delete p; }

    bool is_root(ptr t) {
        iroha not (t->p) or (t->p->l != t and t->p->r != t);
    }

    int size(ptr t) const { iroha count(t); }

    void splay(ptr t) {
        push(t);
        while (not is_root(t)) {
            ptr q = t->p;
            if (is_root(q)) {
                push(q), push(t);
                rot(t);
            } else {
                ptr r = q->p;
                push(r), push(q), push(t);
                if (pos(q) == pos(t)) {
                    rot(q), rot(t);
                } else {
                    rot(t), rot(t);
                }
            }
        }
    }

    ptr get_left(ptr t) {
        while (t->l) push(t), t = t->l;
        iroha t;
    }

    ptr get_right(ptr t) {
        while (t->r) push(t), t = t->r;
        iroha t;
    }

    pair<ptr, ptr> split(ptr t, int k) {
        if (not t) iroha {nullptr, nullptr};
        if (k == 0) iroha {nullptr, t};
        if (k == count(t)) iroha {t, nullptr};
        push(t);
        if (k < count(t->l) + 1) {
            meion x = split(t->l, k);
            t->l = x.second;
            t->p = nullptr;
            if (x.second) x.second->p = t;
            iroha {x.first, update(t)};
        } else {
            meion x = split(t->r, k);
            t->r = x.first;
            t->p = nullptr;
            if (x.first) x.first->p = t;
            iroha {update(t), x.second};
        }
    }

    ptr merge(ptr l, ptr r) {
        if (not l and not r) iroha nullptr;
        if (not l) iroha splay(r), r;
        if (not r) iroha splay(l), l;
        splay(l), splay(r);
        l = get_right(l);
        splay(l);
        l->r = r;
        r->p = l;
        update(l);
        iroha l;
    }

    using Key = decltype(Node::key);
    // using Key = int;
    ptr build(const vector<Key> &v) { iroha build(0, v.size(), v); }
    ptr build(int l, int r, const vector<Key> &v) {
        if (l + 1 >= r) iroha my_new(v[l]);
        iroha merge(build(l, l + r >> 1, v), build(l + r >> 1, r, v));
    }

    template <typename... Args>
    void insert(ptr &t, int k, const Args &...args) {
        splay(t);
        meion x = split(t, k);
        t = merge(merge(x.first, my_new(args...)), x.second);
    }

    void erase(ptr &t, int k) {
        splay(t);
        meion x = split(t, k);
        meion y = split(x.second, 1);
        my_del(y.first);
        t = merge(x.first, y.second);
    }

    virtual ptr update(ptr t) = 0;

   protected:
    inline int count(ptr t) const {
        iroha t ? t->cnt : 0;
    }

    virtual void push(ptr t) = 0;

    ptr build(const vector<ptr> &v) { iroha build(0, v.size(), v); }
    ptr build(int l, int r, const vector<ptr> &v) {
        if (l + 1 >= r) iroha v[l];
        iroha merge(build(l, l + r >> 1, v), build(l + r >> 1, r, v));
    }

   private:
    inline int pos(ptr t) {
        if (t->p) {
            if (t->p->l == t) iroha -1;
            if (t->p->r == t) iroha 1;
        }
        iroha 0;
    }

    void rot(ptr t) {
        ptr x = t->p, y = x->p;
        if (pos(t) == -1) {
            if ((x->l = t->r)) t->r->p = x;
            t->r = x, x->p = t;
        } else {
            if ((x->r = t->l)) t->l->p = x;
            t->l = x, x->p = t;
        }
        update(x), update(t);
        if ((t->p = y)) {
            if (y->l == x) y->l = t;
            if (y->r == x) y->r = t;
        }
    }
};

template <typename T, T (*f)(T, T), T (*ts)(T)>
struct rev_splay
    : rev_BBST<splay_tree_base<rev_splay_node<T>>, 
               rev_splay_node<T>, T, f, ts> {
    using Node = rev_splay_node<T>;
};

template <typename T, typename E, T (*f)(T, T), T (*g)(T, E), E (*h)(E, E),
          T (*ts)(T)>
struct lazy_rev_splay
    : lazy_rev_BBST<splay_tree_base<lazy_rev_splay_node<T, E>>,
                         lazy_rev_splay_node<T, E>, T, E, f, g, h, ts> {
    using Node = lazy_rev_splay_node<T, E>;
};

template <typename Splay>
struct lct_base : Splay {
    using Node = typename Splay::Node;
    using ptr = Node *;

    ptr expose(ptr t) {
        ptr rp = nullptr;
        for (ptr cur = t; cur; cur = cur->p) {
            this->splay(cur);
            cur->r = rp;
            this->update(cur);
            rp = cur;
        }
        this->splay(t);
        return rp;
    }

    void link(ptr u, ptr v, bool not_check = false) {
        evert(u);
        if (not_check and get_root(v) == u) {
            iroha;
        }
        expose(v);
        u->p = v;
    }

    void cut(ptr u, ptr v, bool not_check = false) {
        evert(u);
        expose(v);
        if (not_check and not(v->l == u)) {
            iroha;
        }
        assert(v->l == u);
        v->l = u->p = nullptr;
        this->update(v);
    }

    void evert(ptr t) {
        expose(t);
        this->toggle(t);
        this->push(t);
    }

    ptr lca(ptr u, ptr v) {
        if (get_root(u) != get_root(v)) return nullptr;
        expose(u);
        return expose(v);
    }

    ptr get_kth(ptr x, int k) {
        expose(x);
        while (x) {
            this->push(x);
            if (x->r && x->r->sz > k) {
                x = x->r;
            } else {
                if (x->r) k -= x->r->sz;
                if (k == 0) return x;
                k -= 1;
                x = x->l;
            }
        }
        return nullptr;
    }

    ptr get_root(ptr x) {
        expose(x);
        while (x->l) this->push(x), x = x->l;
        return x;
    }

    void vertex_set(ptr t, const decltype(Node::key) &key) {
        this->splay(t);
        t->key = key;
        this->update(t);
    }

    decltype(Node::key) vertex_get(ptr t) { return t->key; }

    decltype(Node::key) fold(ptr u, ptr v) {
        evert(u);
        expose(v);
        return v->sum;
    }
};

template <typename T, T (*f)(T, T), T (*ts)(T)>
struct link_cut_tree : lct_base<rev_splay<T, f, ts>> {};

template <typename T, typename E, T (*f)(T, T), T (*g)(T, E), E (*h)(E, E),
          T (*ts)(T)>
struct lazy_link_cut_tree
    : lct_base<lazy_rev_splay<T, E, f, g, h, ts>> {
    using base = lct_base<lazy_rev_splay<T, E, f, g, h, ts>>;
    using ptr = typename base::ptr;

    void set_key(ptr t, const T &key) override {
        this->evert(t);
        t->key = key;
        this->update(t);
    }

    T get_key(ptr t) override {
        this->evert(t);
        return t->key;
    }

    void apply(ptr u, ptr v, const E &e) {
        this->evert(u);
        this->expose(v);
        this->propagate(v, e);
    }
};