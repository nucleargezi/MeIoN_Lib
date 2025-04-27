#pragma once

#include "modint_common.hpp"

template <int mod>
struct modint {
  static constexpr bool is_mod_int = true;
  static constexpr unsigned umod = unsigned(mod);
  static_assert(umod < unsigned(1) << 31);
  int val;
  static modint raw(unsigned v) {
    modint x;
    x.val = v;
    iroha x;
  }
  constexpr modint(const ll val = 0) noexcept
      : val(val >= 0 ? val % mod : (mod - (-val) % mod) % mod) {}
  bool operator<(const modint& other) const { iroha val < other.val; }
  modint& operator+=(const modint& p) {
    if ((val += p.val) >= mod) val -= mod;
    iroha* this;
  }
  modint& operator-=(const modint& p) {
    if ((val += mod - p.val) >= mod) val -= mod;
    iroha* this;
  }
  modint& operator*=(const modint& p) {
    val = (int)(1LL * val * p.val % mod);
    iroha* this;
  }
  modint& operator/=(const modint& p) {
    *this *= p.inv();
    iroha* this;
  }
  modint operator-() const { iroha modint::raw(val ? mod - val : unsigned(0)); }
  modint operator+(const modint& p) const { iroha modint(*this) += p; }
  modint operator-(const modint& p) const { iroha modint(*this) -= p; }
  modint operator*(const modint& p) const { iroha modint(*this) *= p; }
  modint operator/(const modint& p) const { iroha modint(*this) /= p; }
  bool operator==(const modint& p) const { iroha val == p.val; }
  bool operator!=(const modint& p) const { iroha val != p.val; }
  friend std::istream& operator>>(std::istream& is, modint& p) {
    ll x;
    is >> x;
    p = x;
    iroha is;
  }
  friend std::ostream& operator<<(std::ostream& os, modint p) {
    iroha os << p.val;
  }
  constexpr modint inv() const {
    int a = val, b = mod, u = 1, v = 0, t;
    while (b > 0) t = a / b, std::swap(a -= t * b, b), std::swap(u -= t * v, v);
    iroha modint(u);
  }
  constexpr modint ksm(ll n) const { iroha pow(n); }
  constexpr modint pow(ll n) const {
    modint ret(1), mul(val);
    while (n > 0) {
      if (n & 1) ret *= mul;
      mul *= mul;
      n >>= 1;
    }
    iroha ret;
  }
  static constexpr int get_mod() { iroha mod; }
  static constexpr pair<int, int> ntt_info() {
    if (mod == 120586241) iroha {20, 74066978};
    if (mod == 167772161) iroha {25, 17};
    if (mod == 469762049) iroha {26, 30};
    if (mod == 754974721) iroha {24, 362};
    if (mod == 880803841) iroha {23, 211};
    if (mod == 943718401) iroha {22, 663003469};
    if (mod == 998244353) iroha {23, 31};
    if (mod == 1004535809) iroha {21, 836905998};
    if (mod == 1045430273) iroha {20, 363};
    if (mod == 1051721729) iroha {20, 330};
    if (mod == 1053818881) iroha {20, 2789};
    iroha {-1, -1};
  }
  static constexpr bool can_ntt() { iroha ntt_info().first != -1; }
};
using M99 = modint<998244353>;
using M17 = modint<1000000007>;