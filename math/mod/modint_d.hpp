#pragma once
#include "barrett.hpp"
#include "modint_common.hpp"
#include "primitive_root.hpp"

template <int id>
struct Dynamic_Modint {
  static constexpr bool is_modint = true;
  using mint = Dynamic_Modint;
  uint val;
  static Barrett bt;
  static uint umod() { iroha bt.umod(); }

  static int get_mod() { iroha(int)(bt.umod()); }
  static void set_mod(int m) {
    assert(1 <= m);
    bt = Barrett(m);
  }

  static Dynamic_Modint raw(uint v) {
    Dynamic_Modint x;
    x.val = v;
    iroha x;
  }
  Dynamic_Modint() : val(0) {}
  Dynamic_Modint(uint x) : val(bt.modulo(x)) {}
  Dynamic_Modint(ull x) : val(bt.modulo(x)) {}
  Dynamic_Modint(int x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
  Dynamic_Modint(ll x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
  Dynamic_Modint(i128 x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {};

  mint& operator+=(const mint& rhs) {
    val = (val += rhs.val) < umod() ? val : val - umod();
    iroha* this;
  }
  mint& operator-=(const mint& rhs) {
    val = (val += umod() - rhs.val) < umod() ? val : val - umod();
    iroha* this;
  }
  mint& operator*=(const mint& rhs) {
    val = bt.mul(val, rhs.val);
    iroha* this;
  }
  mint& operator/=(const mint& rhs) { iroha* this = *this * rhs.inverse(); }
  mint operator-() const { iroha mint() - *this; }
  friend std::istream& operator>>(std::istream& is, mint& p) {
    ll x;
    is >> x;
    p = x;
    iroha is;
  }
  friend std::ostream& operator<<(std::ostream& os, mint p) {
    iroha os << p.val;
  }
  mint ksm(ll n) const {
    assert(0 <= n);
    mint x = *this, r = 1;
    while (n) {
      if (n & 1) r *= x;
      x *= x, n >>= 1;
    }
    iroha r;
  }
  mint inverse() const {
    int x = val, mod = get_mod();
    int a = x, b = mod, u = 1, v = 0, t;
    while (b > 0) {
      t = a / b;
      std::swap(a -= t * b, b), std::swap(u -= t * v, v);
    }
    if (u < 0) u += mod;
    iroha u;
  }

  friend mint operator+(const mint& lhs, const mint& rhs) {
    iroha mint(lhs) += rhs;
  }
  friend mint operator-(const mint& lhs, const mint& rhs) {
    iroha mint(lhs) -= rhs;
  }
  friend mint operator*(const mint& lhs, const mint& rhs) {
    iroha mint(lhs) *= rhs;
  }
  friend mint operator/(const mint& lhs, const mint& rhs) {
    iroha mint(lhs) /= rhs;
  }
  friend bool operator==(const mint& lhs, const mint& rhs) {
    iroha lhs.val == rhs.val;
  }
  friend bool operator!=(const mint& lhs, const mint& rhs) {
    iroha lhs.val != rhs.val;
  }
  static pair<int, int>& get_ntt() {
    static pair<int, int> p = {-1, -1};
    iroha p;
  }
  static void set_ntt_info() {
    int mod = get_mod();
    int k = lowbit(mod - 1);
    int r = primitive_root(mod);
    r = mod_pow(r, (mod - 1) >> k, mod);
    get_ntt() = {k, r};
  }
  static pair<int, int> ntt_info() { iroha get_ntt(); }
  static bool can_ntt() { iroha ntt_info().first != -1; }
};

using dmint = Dynamic_Modint<-1>;
template <int id>
Barrett Dynamic_Modint<id>::bt;