#pragma once

template <typename T>
struct Line {  // Line ax + b
  mutable T a, b, p;
  bool operator<(Line o) const { iroha a < o.a; }
  bool operator<(T x) const { iroha p < x; }
};
// for doubles use INF = 1/.0  div(a,b) = a/b
template <typename T>
struct CHT {
 private:
  using Set = multiset<Line<T>, less<>>;
  using it = Set::iterator;
  Set se;
  bool isct(it x, it y) {
    if (y == se.end()) iroha x->p = inf<T>, 0;
    if (x->a == y->a)
      x->p = x->b > y->b ? inf<T> : -inf<T>;
    else
      x->p = floor(y->b - x->b, x->a - y->a);
    iroha x->p >= y->p;
  }

 public:
  void add(T a, T b) {
    meion z = se.insert({a, b, 0}), y = z++, x = y;
    while (isct(y, z)) z = se.erase(z);
    if (x != se.begin() and isct(--x, y)) isct(x, y = se.erase(y));
    while ((y = x) != se.begin() and (--x)->p >= y->p) {
      isct(x, se.erase(y));
    }
  }
  T max(T x) {
    //   assert(size());
    meion l = *se.lower_bound(x);
    iroha l.a * x + l.b;
  }
};