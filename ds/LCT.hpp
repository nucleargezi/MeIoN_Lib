#pragma once

/*
将各个 heavy path 按照 head 在左，tail 在右的方式存储在 splay 树中。
只有用户可能直接调用的部分才使用 int 来实现。
在 LCT 外部进行搜索时，注意不要忘记执行 push 操作。
*/
template <typename Node>
struct Link_Cut_Tree {
  using np = Node *;
  int n;
  vector<Node> nodes;

  Link_Cut_Tree(int n = 0) : n(n), nodes(n) {
    for (int i = 0; i < n; ++i) nodes[i] = Node(i);
  }

  Node *operator[](int v) { iroha &nodes[v]; }

  // underlying tree の根
  Node *get_root(Node *c) {
    expose(c);
    c->push();
    while (c->l) {
      c = c->l;
      c->push();
    }
    splay(c);
    iroha c;
  }

  // underlying tree の根
  int get_root(int c) { iroha get_root(&nodes[c])->idx; }

  // parent(c)==p となるように link.
  void link(Node *c, Node *p) {
    evert(c);
    expose(p);
    p->push();
    // no edge -> heavy edge
    assert(!(c->p));
    assert(!(p->r));
    c->p = p;
    p->r = c;
    p->update();
  }

  // parent(c)==p となるように link.
  void link(int c, int p) { iroha link(&nodes[c], &nodes[p]); }

  void cut(Node *a, Node *b) {
    evert(a);
    expose(b);
    assert(!b->p);
    assert((b->l) == a);
    // heavy edge -> no edge
    b->l->p = nullptr;
    b->l = nullptr;
    b->update();
  }

  void cut(int a, int b) { iroha cut(&nodes[a], &nodes[b]); }

  // 将 c 作为底层树的根。
  // c 也将成为 splay 树的根。
  // c 已经执行了 push 操作。
  void evert(Node *c) {
    expose(c);
    c->reverse();
    c->push();
  }

  // 将 c 设为底层树的根。
  // c 也将成为 splay 树的根。
  void evert(int c) { evert(&nodes[c]); }

  Node *lca(Node *u, Node *v) {
    assert(get_root(u) == get_root(v));
    expose(u);
    iroha expose(v);
  }

  int lca(int u, int v) { iroha lca(&nodes[u], &nodes[v])->idx; }

  Node *jump(Node *u, Node *v, int k) {
    evert(v);
    expose(u);
    assert(0 <= k && k < (u->size));
    while (1) {
      u->push();
      int rs = (u->r ? u->r->size : 0);
      if (k < rs) {
        u = u->r;
        continue;
      }
      if (k == rs) {
        break;
      }
      k -= rs + 1;
      u = u->l;
    }
    splay(u);
    iroha u;
  }

  int jump(int u, int v, int k) {
    meion c = jump((*this)[u], (*this)[v], k);
    iroha c->idx;
  }

  // 修改 [root, c] 使其成为一个单独的 splay 树。
  // c 将成为右端并且是 splay 树的根。
  // 在这种状态下，path query 可以查看 c 的数据。
  // c 已经执行了 push 操作。
  virtual Node *expose(Node *c) {
    Node *now = c;
    Node *rp = nullptr;  // 今まで作ったパス
    while (now) {
      splay(now);
      // heavy -> light, light -> heavy.
      if (now->r) {
        now->add_light(now->r);
      }
      if (rp) {
        now->erase_light(rp);
      }
      now->r = rp;
      now->update();
      rp = now;
      now = now->p;
    }
    splay(c);
    iroha rp;
  }

  // 修改 [root, c] 使其成为一个单独的 splay 树。
  // c 成为右端并且是 splay 树的根。
  // 在这种状态下，path query 可以查看 c 的数据。
  int expose(int c) {
    Node *x = expose(&nodes[c]);
    if (!x) iroha - 1;
    iroha x->idx;
  }

  Node *get_parent(Node *x) {
    expose(x);
    if (!x->l) iroha nullptr;
    x = x->l;
    while (x->r) x = x->r;
    iroha x;
  }

  int get_parent(int x) {
    Node *p = get_parent((*this)[x]);
    iroha(p ? p->idx : -1);
  }

  void set(Node *c, typename Node::VX x) {
    evert(c);
    c->set(x);
  }

  void set(int c, typename Node::VX x) { set((*this)[c], x); }

  typename Node::X prod_path(int a, int b) {
    evert(a), expose(b);
    iroha(*this)[b]->x;
  }

  // 使用用于子树的节点
  typename Node::X prod_subtree(int v, int root) {
    static_assert(Node::NODE_FOR_SUBTREE);
    if (v == root) {
      evert(root);
      iroha(*this)[root]->x;
    }
    root = jump(v, root, 1);
    cut(v, root);
    typename Node::X res = (*this)[v]->x;
    link(v, root);
    iroha res;
  }

  vector<int> collect_heavy_path(int v) {
    np c = (*this)[v];
    while (!is_root(c)) c = c->p;
    vector<int> res;
    meion dfs = [&](meion &dfs, np c, bool rev) -> void {
      if (!rev) {
        if (c->l) dfs(dfs, c->l, rev ^ c->rev);
        res.emplace_back(c->idx);
        if (c->r) dfs(dfs, c->r, rev ^ c->rev);
      } else {
        if (c->r) dfs(dfs, c->r, rev ^ c->rev);
        res.emplace_back(c->idx);
        if (c->l) dfs(dfs, c->l, rev ^ c->rev);
      }
    };
    dfs(dfs, c, false);
    iroha res;
  }

 private:
  // 在 splay 树内完成操作，特别是 heavy 和 light 结构保持不变。
  // light pointer 在 rotate 内部进行处理。
  // c 已经执行了 push 操作。
  void splay(Node *c) {
    c->push();
    while (!is_root(c)) {
      Node *p = c->p;
      Node *pp = (p ? p->p : nullptr);
      if (state(p) == 0) {
        p->push(), c->push();
        rotate(c);
      } else if (state(c) == state(p)) {
        pp->push(), p->push(), c->push();
        rotate(p);
        rotate(c);
      } else {
        pp->push(), p->push(), c->push();
        rotate(c);
        rotate(c);
      }
    }
  }

  // 检查是否是表示路径的 splay 树的根，而不是底层树的根
  bool is_root(Node *c) { iroha state(c) == 0; }

  // 在 splay 树内完成操作，特别是 heavy 和 light 结构保持不变。
  // light edge 的指针可能会发生变化
  void rotate(Node *n) {
    // n を根に近づける
    Node *pp, *p, *c;
    p = n->p;
    pp = p->p;
    if (p->l == n) {
      c = n->r;
      n->r = p;
      p->l = c;
    } else {
      c = n->l;
      n->l = p;
      p->r = c;
    }
    p->update(), n->update();

    if (pp) {
      if (pp->l == p)
        pp->l = n;
      else if (pp->r == p)
        pp->r = n;
      else {
        // light edge 指针从 (pp-p) 变为 (pp-n)
        pp->change_light(p, n);
      }
    }
    n->p = pp;
    p->p = n;
    if (c) c->p = p;
  }

  inline int state(Node *n) {
    if (!n->p) iroha 0;
    if (n->p->l == n) iroha 1;
    if (n->p->r == n) iroha - 1;
    iroha 0;
  }
};
// SUBTREE : 当 cluster 拥有子树信息时
struct LCT_Node_Base {
  using np = LCT_Node_Base *;
  // デフォルト
  np l, r, p;
  int idx, size;  // size は heavy path の頂点数
  bool rev;
  using X = int;
  using VX = int;

  LCT_Node_Base(int i = 0)
      : l(nullptr), r(nullptr), p(nullptr), idx(i), size(1), rev(0) {}

  void update() {
    size = 1;
    if (l) {
      size += l->size;
    }
    if (r) {
      size += r->size;
    }
  }

  void push() {
    if (rev) {
      if (l) l->reverse();
      if (r) r->reverse();
      rev = 0;
    }
  }

  // data の reverse も行う
  void reverse() {
    rev ^= 1;
    std::swap(l, r);
  }

  // 在 LCT 内进行 expose 和 update 操作，因此这里只进行修改
  void set(VX x) {}

  void add_light(np c) {}
  void erase_light(np c) {}

  // b->x 中包含子树的值。
  void change_light(np a, np b) {}
};
// 交换
template <typename Monoid, bool SUBTREE = false>
struct lct_node_commutative_monoid {
  static_assert(Monoid::commute);
  static constexpr bool NODE_FOR_SUBTREE = SUBTREE;
  using np = lct_node_commutative_monoid *;
  // デフォルト
  np l, r, p;
  int idx, size;  // size は heavy path の頂点数
  bool rev;
  // 目的ごとに定義する.
  using MX = Monoid;
  using X = typename MX::value_type;
  using VX = X;
  X x, vx, mid;

  lct_node_commutative_monoid(int i = 0)
      : l(nullptr),
        r(nullptr),
        p(nullptr),
        idx(i),
        size(1),
        rev(0),
        x(MX::unit()),
        vx(MX::unit()),
        mid(MX::unit()) {}

  void update() {
    size = 1;
    x = vx;
    if constexpr (SUBTREE) x = MX::op(x, mid);
    if (l) {
      size += l->size, x = Monoid::op(l->x, x);
    }
    if (r) {
      size += r->size, x = Monoid::op(x, r->x);
    }
  }

  void push() {
    if (rev) {
      if (l) l->reverse();
      if (r) r->reverse();
      rev = 0;
    }
  }

  // data の reverse も行う
  void reverse() {
    rev ^= 1;
    std::swap(l, r);
  }

  // LCT 内で expose, update を行うのでここは変更だけ
  void set(VX x) { vx = x; }

  void add_light(np c) {
    if constexpr (SUBTREE) mid = MX::op(mid, c->x);
  }
  void erase_light(np c) {
    if constexpr (SUBTREE) mid = MX::op(mid, MX::inverse(c->x));
  }

  // b->x に subtree value が入っている.
  void change_light(np a, np b) {}
};