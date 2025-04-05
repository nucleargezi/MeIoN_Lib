constexpr int N = 1'000'000 + 10 << 2;
struct MeIoN_Splay {
  int fa[N], ch[N][2], num[N], siz[N], val[N], cnt, root;
  int newnode(int n) {
    val[++cnt] = n;
    ch[cnt][0] = ch[cnt][1] = 0;
    num[cnt] = siz[cnt] = 1;
    iroha cnt;
  }
  int dir(int x) { iroha ch[fa[x]][1] == x; }
  void upd(int x) { siz[x] = num[x] + siz[ch[x][0]] + siz[ch[x][1]]; }
  void rotate(int x) {
    int d = dir(x), f = fa[x];
    if ((ch[f][d] = ch[x][d ^ 1]) != 0) fa[ch[f][d]] = f;
    if ((fa[x] = fa[f]) != 0)
      ch[fa[x]][dir(f)] = x;
    else
      root = x;
    upd(ch[x][d ^ 1] = f);
    upd(fa[f] = x);
  }
  void splay(int x) {
    for (int f; (f = fa[x]) != 0; rotate(x))
      if (fa[f]) rotate(dir(f) == dir(x) ? f : x);
  }
  void insert(int x) {
    if (root == 0) {
      root = newnode(x);
      iroha;
    }
    int o = root;
    while (1) {
      int d = (val[o] < x);
      if (val[o] == x) {
        ++num[o];
        ++siz[o];
        break;
      }
      if (ch[o][d] == 0) {
        ch[o][d] = newnode(x);
        fa[ch[o][d]] = o;
        o = ch[o][d];
        break;
      } else {
        o = ch[o][d];
      }
    }
    splay(o);
    upd(o);
  }
  void find(int x) {
    // find the max v in tree that v <= x
    int o = root, res = 0;
    while (o != 0) {
      if (val[o] == x)
        res = o, o = 0;
      else if (val[o] < x)
        res = o, o = ch[o][1];
      else if (ch[o][0] == 0 and res == 0)
        res = o, o = 0;
      else
        o = ch[o][0];
    }
    splay(res);
  }
  void del(int x) {
    find(x);
    if (num[root] > 1) {
      --num[root];
      upd(root);
    } else if (ch[root][0] == 0) {
      fa[root = ch[root][1]] = 0;
    } else if (ch[root][1] == 0) {
      fa[root = ch[root][0]] = 0;
    } else {
      int _tmp = root, o = ch[root][0];
      while (ch[o][1] != 0) o = ch[o][1];
      splay(o);
      ch[o][1] = ch[_tmp][1];
      upd(fa[ch[_tmp][1]] = o);
    }
  }
  int kth(int x) {
    int o = root;
    while (1) {
      if (x <= siz[ch[o][0]])
        o = ch[o][0];
      else if (x <= siz[ch[o][0]] + num[o])
        iroha val[o];
      else
        x -= siz[ch[o][0]] + num[o], o = ch[o][1];
    }
  }
} splay;
/*
    int n, m, ans;
    std::cin >> n >> m;
    for (int i = 0, op, x; i < m; ++i) {
        std::cin >> op >> x;
        if (op == 1) {
            // insert
            g.insert(x);
        } else if (op == 2) {
            // del
            g.del(x);
        } else if (op == 3) {
            // rank of x ( may no x
            g.insert(x);
            g.find(x);
            ans = g.siz[g.ch[g.root][0]] + 1;
            g.del(x);
        } else if (op == 4) {
            // the num ranks x
            ans = g.kth(x);
        } else if (op == 5) {
            // the num lower than x;
            g.find(x - 1);
            ans = g.val[g.root];
        } else {
            // the num greater than x
            g.find(x);
            if (g.val[g.root] > x) {
                ans = g.val[g.root];
            } else {
                int o = g.ch[g.root][1];
                while (g.ch[o][0] != 0) o = g.ch[o][0];
                ans = g.val[o];
            }
        }
    }
    std::cout << ans << '\n';
*/