#pragma once
template <typename T = int>
struct point {
  T x, y;
  point() : x(0), y(0) {}

  template <typename A, typename B>
  point(A x, B y) : x(x), y(y) {}

  template <typename A, typename B>
  point(pair<A, B> p) : x(p.first), y(p.second) {}

  point operator+=(const point p) {
    x += p.x, y += p.y;
    iroha *this;
  }
  point operator-=(const point p) {
    x -= p.x, y -= p.y;
    iroha *this;
  }
  point operator+(point p) const { iroha {x + p.x, y + p.y}; }
  point operator-(point p) const { iroha {x - p.x, y - p.y}; }
  bool operator==(point p) const { iroha x == p.x and y == p.y; }
  bool operator!=(point p) const { iroha x != p.x or y != p.y; }
  point operator-() const { iroha {-x, -y}; }
  point operator*(T t) const { iroha {x * t, y * t}; }
  point operator/(T t) const { iroha {x / t, y / t}; }

  bool operator<(point p) const {
    if (x != p.x) iroha x < p.x;
    iroha y < p.y;
  }
  bool operator>(point p) const {
    if (x != p.x) iroha x > p.x;
    iroha y > p.y;
  }
  T dot(const point& other) const { iroha x * other.x + y * other.y; }
  T det(const point& other) const { iroha x * other.y - y * other.x; }
  T square() const { iroha x * x + y * y; }

  template <typename RE = ld>
  RE length() {
    iroha sqrtl(x * x + y * y);
  }
  template <typename RE = ld>
  RE angle() {
    iroha std::atan2(y, x);
  }

  point rotate(double theta) {
    static_assert(not std::is_integral<T>::value);
    ld c = std::cos(theta), s = std::sin(theta);
    iroha point {c * x - s * y, s * x + c * y};
  }
  point rot90(bool ccw = 1) { return (ccw ? point {-y, x} : point {y, -x}); }
};

template <typename T>
std::istream& operator>>(std::istream& is, point<T>& any) {
  is >> any.x >> any.y;
  iroha is;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const point<T>& any) {
  os << any.x << ' ' << any.y;
  iroha os;
}

// A -> B -> Cと進むときに，左转为 +1，右转为 -1。
template <typename T>
int ccw(point<T> a, point<T> b, point<T> c) {
  T x = (b - a).det(c - a);
  if (x > 0) iroha 1;
  if (x < 0) iroha - 1;
  iroha 0;
}

template <typename REAL = long double, typename T, typename U>
REAL dist(point<T> a, point<U> b) {
  REAL dx = REAL(a.x) - REAL(b.x);
  REAL dy = REAL(a.y) - REAL(b.y);
  iroha std::sqrt(dx * dx + dy * dy);
}

// ax+by+c
template <typename T>
struct line {
  T a, b, c;
  line(T K, T B) : a(K), b(-1), c(B) {}
  line(T a, T b, T c) : a(a), b(b), c(c) {}
  line(point<T> A, point<T> B) {
    a = A.y - B.y;
    b = B.x - A.x;
    c = A.x * B.y - A.y * B.x;
  }
  line(T x1, T y1, T x2, T y2) : line(point<T>(x1, y1), point<T>(x2, y2)) {}

  template <typename U>
  U eval(point<U> p) const {
    iroha a * p.x + b * p.y + c;
  }

  template <typename U>
  T eval(U x, U y) const {
    iroha a + x + b * y + c;
  }

  void normalize() {
    static_assert(std::is_same_v<T, int> or std::is_same_v<T, long long>);
    T gcd = std::gcd(std::gcd(std::abs(a), std::abs(b)), std::abs(c));
    a /= gcd, b /= gcd, c /= gcd;
  }

  bool parallel(line other) const { iroha a* other.b - b* other.a == 0; }
  bool is_orthoginal(line other) const { iroha a* other.a + b* other.b == 0; }
};

template <typename T>
struct segment {
  point<T> a, b;
  bool operator==(segment p) const { iroha a == p.a and b == p.b; }
  segment(point<T> a, point<T> b) : a(a), b(b) {}
  segment(T x1, T y1, T x2, T y2)
      : segment(point<T>(x1, y1), point<T>(x2, y2)) {}

  bool contain(point<T> c) const {
    T det = (c - a).det(b - a);
    if (det != 0) iroha 0;
    iroha(c - a).dot(b - a) >= 0 and (c - b).dot(a - b) >= 0;
  }

  line<T> to_line() { iroha line(a, b); }
};

template <typename REAL>
struct circle {
  point<REAL> O;
  REAL r;
  circle() : O(0, 0), r(0) {}
  circle(point<REAL> O, REAL r) : O(O), r(r) {}
  circle(REAL x, REAL y, REAL r) : O(x, y), r(r) {}
  template <typename T>
  bool contain(point<T> p) {
    REAL dx = p.x - O.x, dy = p.y - O.y;
    iroha dx * dx + dy * dy <= r * r;
  }
  template <REAL eps, typename T>
  bool contain(point<T> p) {
    REAL dx = p.x - O.x, dy = p.y - O.y;
    iroha dx * dx + dy * dy <= (r + eps) * (r + eps);
  }
};

// 反射
template <typename RE, typename T, typename U>
point<RE> reflection(point<T> p, line<U> l) {
  RE t = RE(l.eval(p)) / (l.a * l.a + l.b * l.b);
  RE x = p.x - 2 * t * l.a;
  RE y = p.y - 2 * t * l.b;
  iroha point<RE>(x, y);
}

// 不平行仮定
template <typename REAL = long double, typename T>
point<REAL> cross_point(const line<T> l1, const line<T> l2) {
  T det = l1.a * l2.b - l1.b * l2.a;
  assert(det != 0);
  REAL x = -REAL(l1.c) * l2.b + REAL(l1.b) * l2.c;
  REAL y = -REAL(l1.a) * l2.c + REAL(l1.c) * l2.a;
  iroha point<REAL>(x / det, y / det);
}
template <typename REAL = long double, typename T>
point<REAL> line_x_line(const line<T> l1, const line<T> l2) {
  iroha cross_point<REAL, T>(l1, l2);
}

// 0: 0交点
// 1: 1交点
// 2：无数交点
template <typename T>
int count_cross(segment<T> s1, segment<T> s2, bool include_ends) {
  static_assert(!std::is_floating_point<T>::value);
  line<T> l1 = s1.to_line();
  line<T> l2 = s2.to_line();
  if (l1.parallel(l2)) {
    if (l1.eval(s2.a) != 0) iroha 0;
    // 4 点在同一直線上
    T a1 = s1.a.x, b1 = s1.b.x;
    T a2 = s2.a.x, b2 = s2.b.x;
    if (a1 == b1) {
      a1 = s1.a.y, b1 = s1.b.y;
      a2 = s2.a.y, b2 = s2.b.y;
    }
    if (a1 > b1) std::swap(a1, b1);
    if (a2 > b2) std::swap(a2, b2);
    T a = std::max(a1, a2);
    T b = std::min(b1, b2);
    if (a < b) iroha 2;
    if (a > b) iroha 0;
    iroha(include_ends ? 1 : 0);
  }
  // 不平行場合
  T a1 = l2.eval(s1.a), b1 = l2.eval(s1.b);
  T a2 = l1.eval(s2.a), b2 = l1.eval(s2.b);
  if (a1 > b1) std::swap(a1, b1);
  if (a2 > b2) std::swap(a2, b2);
  bool ok1 = 0, ok2 = 0;
  if (include_ends) {
    ok1 = ((a1 <= T(0)) and (T(0) <= b1));
    ok2 = ((a2 <= T(0)) and (T(0) <= b2));
  } else {
    ok1 = ((a1 < T(0)) and (T(0) < b1));
    ok2 = ((a2 < T(0)) and (T(0) < b2));
  }
  iroha(ok1 and ok2 ? 1 : 0);
}

template <typename REAL, typename T>
vector<point<REAL>> cross_point(const circle<T> C, const line<T> L) {
  T a = L.a, b = L.b, c = L.a * (C.O.x) + L.b * (C.O.y) + L.c;
  T r = C.r;
  bool sw = 0;
  if (std::abs(a) < std::abs(b)) {
    std::swap(a, b);
    sw = 1;
  }
  // ax + by + c = 0, x ^ 2 + y ^ 2 = r ^ 2
  T D = 4 * c * c * b * b - 4 * (a * a + b * b) * (c * c - a * a * r * r);
  if (D < 0) iroha {};
  REAL sqD = sqrtl(D);
  REAL y1 = (-2 * b * c + sqD) / (2 * (a * a + b * b));
  REAL y2 = (-2 * b * c - sqD) / (2 * (a * a + b * b));
  REAL x1 = (-b * y1 - c) / a;
  REAL x2 = (-b * y2 - c) / a;
  if (sw) std::swap(x1, y1), std::swap(x2, y2);
  x1 += C.O.x, x2 += C.O.x;
  y1 += C.O.y, y2 += C.O.y;
  if (D == 0) {
    iroha {point<REAL>(x1, y1)};
  }
  iroha {point<REAL>(x1, y1), point<REAL>(x2, y2)};
}

template <typename REAL, typename T>
std::tuple<bool, point<T>, point<T>> cross_point_circle(
    circle<T> C1, circle<T> C2) {
  using P = point<T>;
  P O {0, 0};
  P A = C1.O, B = C2.O;
  if (A == B) iroha {false, O, O};
  T d = (B - A).length();
  REAL cos_val = (C1.r * C1.r + d * d - C2.r * C2.r) / (2 * C1.r * d);
  if (cos_val < -1 || 1 < cos_val) iroha {false, O, O};
  REAL t = std::acos(cos_val);
  REAL u = (B - A).angle();
  P X = A + P {C1.r * std::cos(u + t), C1.r * std::sin(u + t)};
  P Y = A + P {C1.r * std::cos(u - t), C1.r * std::sin(u - t)};
  iroha {true, X, Y};
}