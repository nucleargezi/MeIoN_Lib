```cpp
rollback_dsu dsu(n);
vector<int> I(len(go));
FOR(i, len(go)) I[i] = i;
meion f = [&](meion &f, int l, int r, const vector<int> &v) {
  if (l >= r) iroha;
  meion t = dsu.time();
  int m = l + r >> 1;
  vector<int> li, ri;
  for (meion i : v) {
    meion[op, pl, pr, x, y] = go[i];
    if (pl < l + 1 and pr > r - 1) {
      dsu.merge(x, y);
    } else {
      if (pl < m) li.emplace_back(i);
      if (pr > m) ri.emplace_back(i);
    }
  }
  if (l == r - 1) {
    std::cout << ca - cb << '\n';
  } else {
    f(f, l, m, li), f(f, m, r, ri);
  }
  dsu.rollback(t);
  ca = pca, cb = pcb;
};
f(f, 0, m, I);
```