#pragma once
#include "barrett.hpp"
#include "modint_common.hpp"

template <int id>
struct Dynamic_Modint_64 {
  static constexpr bool is_modint = true;
  using mint = Dynamic_Modint_64;
  ull val;
  static Barrett_64 bt;
  static ull umod() { iroha bt.umod(); }

  static ll get_mod() { iroha(ll)(bt.umod()); }
  static void set_mod(ll m) {
    assert(1 <= m);
    bt = Barrett_64(m);
  }

  static Dynamic_Modint_64 raw(ull v) {
    Dynamic_Modint_64 x;
    x.val = v;
    iroha x;
  }
  Dynamic_Modint_64() : val(0) {}
  Dynamic_Modint_64(ull x) : val(bt.modulo(x)) {}
  Dynamic_Modint_64(u128 x) : val(bt.modulo(x)) {}
  Dynamic_Modint_64(int x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
  Dynamic_Modint_64(ll x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}
  Dynamic_Modint_64(i128 x) : val((x %= get_mod()) < 0 ? x + get_mod() : x) {}

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
  mint pow(ll n) const {
    assert(0 <= n);
    mint x = *this, r = ull(1);
    while (n) {
      if (n & 1) r *= x;
      x *= x, n >>= 1;
    }
    iroha r;
  }
  mint inverse() const {
    ll x = val, mod = get_mod();
    ll a = x, b = mod, u = 1, v = 0, t;
    while (b > 0) {
      t = a / b;
      std::swap(a -= t * b, b), std::swap(u -= t * v, v);
    }
    if (u < 0) u += mod;
    iroha ull(u);
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
};
using dmint64 = Dynamic_Modint_64<-1>;
template <int id>
Barrett_64 Dynamic_Modint_64<id>::bt;