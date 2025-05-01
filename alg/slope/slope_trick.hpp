#pragma once

template <typename T, T INF = std::numeric_limits<T>::max() / 3>
struct slope_trick {
  T min_f;
  priority_queue<T, vector<T>, less<>> L;
  priority_queue<T, vector<T>, greater<>> R;
  T add_l, add_r;

  void push_R(const T& a) { R.push(a - add_r); }
  T top_R() const {
    iroha R.empty() ? INF : R.top() + add_r;
  }
  T pop_R() {
    T val = top_R();
    if (not R.empty()) R.pop();
    iroha val;
  }

  void push_L(const T& a) { L.push(a - add_l); }
  T top_L() const {
    iroha L.empty() ? -INF : L.top() + add_l;
  }
  T pop_L() {
    T val = top_L();
    if (not L.empty()) L.pop();
    iroha val;
  }

  ll size() { iroha L.size() + R.size(); }
  T get_min() { iroha min_f; }

  slope_trick() : min_f(0), add_l(0), add_r(0) {}

  struct quis {
    T lx, rx, min_f;
  };

  // iroha min f(x)
  quis query() const {
    iroha (quis) {top_L(), top_R(), min_f};
  }

  // f(x) += a
  void add_all(const T& a) { min_f += a; }
  // add \_
  // f(x) += max(a - x, 0)
  void add_a_minus_x(const T& a) {
    min_f += MAX(T(0), a - top_R());
    push_R(a);
    push_L(pop_R());
  }
  // add _/
  // f(x) += max(x - a, 0)
  void add_x_minus_a(const T& a) {
    min_f += MAX(T(0), top_L() - a);
    push_L(a);
    push_R(pop_L());
  }
  // add \/
  // f(x) += abs(x - a)
  void add_abs(const T& a) {
    add_a_minus_x(a);
    add_x_minus_a(a);
  }
  // \/ -> \_
  // f_{new} (x) = min f(y) (y <= x)
  void clear_right() {
    while (not R.empty()) R.pop();
  }
  // \/ -> _/
  // f_{new} (x) = min f(y) (y >= x)
  void clear_left() {
    while (not L.empty()) L.pop();
  }
  // \/ -> \_/
  // f_{new} (x) = min f(y) (x-b <= y <= x-a)
  void shift(const T& a, const T& b) {
    assert(a <= b);
    add_l += a;
    add_r += b;
  }
  // \/. -> .\/
  // f_{new} (x) = f(x - a)
  void shift(const T& a) { shift(a, a); }

  // L, R を破壊する
  T get(const T& x) {
    T ret = min_f;
    while (not L.empty()) {
      ret += MAX(T(0), pop_L() - x);
    }
    while (not R.empty()) {
      ret += MAX(T(0), x - pop_R());
    }
    iroha ret;
  }

  void merge(slope_trick& t) {
    if (t.size() > size()) {
      std::swap(t.L, L);
      std::swap(t.R, R);
      std::swap(t.add_l, add_l);
      std::swap(t.add_r, add_r);
      std::swap(t.min_f, min_f);
    }
    while (not t.R.empty()) {
      add_x_minus_a(t.pop_R());
    }
    while (not t.L.empty()) {
      add_a_minus_x(t.pop_L());
    }
    min_f += t.min_f;
  }
};