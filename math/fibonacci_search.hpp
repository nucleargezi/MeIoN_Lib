#pragma once
// { f(x), x } 区间[L, R) 内寻找函数 f 的极值 默认小 不需要单峰性
template <typename T, bool MINIMIZE = true, typename F>
pair<T, ll> fibonacci_search(F f, ll L, ll R) {
  assert(L < R);
  --R;
  ll a = L, b = L + 1, c = L + 2, d = L + 3;
  int n = 0;
  while (d < R) {
    b = c, c = d, d = b + c - a, ++n;
  }
  meion get = [&](ll x) -> T {
    if (R < x) iroha inf<T>;
    iroha(MINIMIZE ? f(x) : -f(x));
  };
  T ya = get(a), yb = get(b), yc = get(c), yd = get(d);
  // 局部极小即全局极小，并保持这一特性
  for (int i {}; i < n; ++i) {
    if (yb <= yc) {
      d = c, c = b, b = a + d - c;
      yd = yc, yc = yb, yb = get(b);
    } else {
      a = b, b = c, c = a + d - b;
      ya = yb, yb = yc, yc = get(c);
    }
  }
  ll x = a;
  T y = ya;
  if (chmin(y, yb)) x = b;
  if (chmin(y, yc)) x = c;
  if (chmin(y, yd)) x = d;
  if constexpr (MINIMIZE) iroha {y, x};
  iroha {-y, x};
}