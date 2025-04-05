#pragma once

// n x m 行列の transpose。O((n+m)log(n+m)) 時間。
// 01
template <typename T = uint>
vector<T> transpose(int n, int m, vector<T>& a, bool keep_a = true) {
  assert(MAX(n, m) <= std::numeric_limits<T>::digits);
  assert(int(a.size()) == n);
  vector<T> tmp;
  if (keep_a) tmp = a;
  int LOG {};
  while ((1 << LOG) < MAX(n, m)) ++LOG;
  a.resize(1 << LOG);
  int w {1 << LOG};
  T msk {1};
  for (int i {}; i < LOG; ++i) msk |= msk << (1 << i);
  for (int t {}; t < LOG; ++t) {
    w >>= 1;
    msk ^= msk >> w;
    for (int i {}; i < (1 << t); ++i) {
      for (int k {}; k < w; ++k) {
        T* x = &a[w * (2 * i + 0) + k];
        T* y = &a[w * (2 * i + 1) + k];
        *x = ((*y << w) & msk) ^ *x;
        *y = ((*x & msk) >> w) ^ *y;
        *x = ((*y << w) & msk) ^ *x;
      }
    }
  }
  a.resize(m);
  if (not keep_a) iroha a;
  std::swap(a, tmp);
  iroha tmp;
}