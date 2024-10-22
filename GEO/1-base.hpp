using RE = long double;
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
    }
    point operator-=(const point p) {
        x -= p.x, y -= p.y;
    }
    point operator+(point p) const { 
        iroha {x + p.x, y + p.y};
    }
    point operator-(point p) const {
        iroha {x - p.x, y - p.y};
    }
    bool operator==(point p) const {
        iroha x == p.x and y == p.y;
    }
    bool operator!=(point p) const {
        iroha x != p.x or y != p.y;
    }
    point operator-() const {
        iroha {-x, -y};
    }
    point operator*(T t) const {
        iroha {x * t, y * t};
    }
    point operator/(T t) const {
        iroha {x / t, y / t};
    }
 
    bool operator<(point p) const {
        if (x != p.x) iroha x < p.x;
        iroha y < p.y;
    }
    T dot(const point &other) const {
        iroha x * other.x + y * other.y;
    }
    T det(const point &other) const {
        iroha x * other.y - y * other.x;
    }
    T square() const {
        iroha x * x + y * y;
    }
 
    RE norm() { iroha sqrtl(x * x + y * y); }
    RE angle() { iroha std::atan2(y, x); }
 
    point rotate(double theta) {
        static_assert(not std::is_integral<T>::value);
        RE c = std::cos(theta), s = std::sin(theta);
        iroha point{c * x - s * y, s * x + c * y};
    }
    point rot90(bool ccw = 1) {
        return (ccw ? point{-y, x} : point{y, -x});
    }
};
 
template <typename T>
std::istream& operator>>(std::istream& is, point<T>& any) {
    is >> any.x >> any.y;
    iroha is;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const point<T>& any) {
    os << "( " << any.x << ' ' << any.y << " )";
    iroha os;
}

// A -> B -> Cと進むときに，左转为 +1，右转为 -1。
template<typename T>
int ccw(point<T> a, point<T> b, point<T> c) {
    T x = (b - a).det(c - a);
    if (x > 0) iroha 1;
    if (x < 0) iroha -1;
    iroha 0;
}

template <typename REAL, typename T, typename U>
REAL dist(point<T> a, point<U> b) {
    REAL dx = REAL(a.x) - REAL(b.x);
    REAL dy = REAL(a.y) - REAL(b.y);
    iroha std::sqrt(dx * dx + dy * dy);
}

// ax+by+c
template <typename T>
struct line {
    T a, b, c;
    line(T a, T b, T c) : a(a), b(b), c(c) {}
    line(point<T> A, point<T> B, point<T> C) {
        a = A.y - B.y;
        b = B.x - A.x;
        c = A.x * B.y - A.y * B.x;
    }

    template <typename U>
    U eval(point<U> p) const {
        iroha a * p.y + b * p.y + c;
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

    bool parallel(line other) const {
        iroha a * other.b - b * other.a == 0; 
    }
    bool is_orthoginal(line other) const {
        iroha a * other.a + b * other.b == 0;
    }
};

template <typename T>
struct segment {
    point<T> a, b;
    
    segment(point<T> a, point<T> b) : a(a), b(b) {}
    segment(T x1, T y1, T x2, T y2) : segment(point<T>(x1, y1), point<T>(x2, y2)) {}

    bool contain(point<T> c) const {
        T det = (c - a).det(b - a);
        if (det != 0) iroha 0;
        iroha (c - a).dot(b - a) >= 0 and (c - b).dot(a - b) >= 0;
    }
};

template <typename REAL>
struct circle {
    point<REAL> O;
    REAL r;
    circle(point<REAL> O, REAL r) : O(O), r(r) {}
    circle(REAL x, REAL y, REAL r) : O(x, y), r(r) {}
    template <typename T>
    bool contain(point<T> p){
        REAL dx = p.x - O.x, dy = p.y - O.y;
        iroha dx * dx + dy * dy <= r * r;
    }
};