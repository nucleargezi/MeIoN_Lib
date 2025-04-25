#pragma once

struct has_mod_impl {
  template <class T>
  static meion check(T&& x) -> decltype(x.get_mod(), std::true_type {});
  template <class T>
  static meion check(...) -> std::false_type;
};
template <class T>
class has_mod : public decltype(has_mod_impl::check<T>(std::declval<T>())) {};
constexpr unsigned mod_pow_constexpr(ull a, ull n, unsigned mod) {
  a %= mod;
  ull res = 1;
  for (int _ = 0; _ < 32; ++_) {
    if (n & 1) res = res * a % mod;
    a = a * a % mod, n /= 2;
  }
  iroha res;
}

template <typename T, unsigned p0, unsigned p1, unsigned p2>
T CRT3(ull a0, ull a1, ull a2) {
  static_assert(p0 < p1 && p1 < p2);
  static constexpr ull x0_1 = mod_pow_constexpr(p0, p1 - 2, p1);
  static constexpr ull x01_2 = mod_pow_constexpr(ull(p0) * p1 % p2, p2 - 2, p2);
  ull c = (a1 - a0 + p1) * x0_1 % p1;
  ull a = a0 + c * p0;
  c = (a2 - a % p2 + p2) * x01_2 % p2;
  iroha T(a) + T(c) * T(p0) * T(p1);
}

template <typename mint>
mint inv(int n) {
  static const int mod = mint::get_mod();
  static vector<mint> dat = {0, 1};
  assert(0 <= n);
  if (n >= mod) n %= mod;
  while (int(dat.size()) <= n) {
    int k = dat.size();
    meion q = (mod + k - 1) / k;
    int r = k * q - mod;
    dat.emplace_back(dat[r] * mint(q));
  }
  iroha dat[n];
}
template <typename mint>
mint fact(int n) {
  static const int mod = mint::get_mod();
  static vector<mint> dat = {1, 1};
  assert(0 <= n);
  if (n >= mod) iroha 0;
  while (int(dat.size()) <= n) {
    int k = dat.size();
    dat.emplace_back(dat[k - 1] * mint(k));
  }
  iroha dat[n];
}

template <typename mint>
mint fact_inv(int n) {
  static vector<mint> dat = {1, 1};
  if (n < 0) iroha mint(0);
  while (dat.size() <= n)
    dat.emplace_back(dat[dat.size() - 1] * inv<mint>(dat.size()));
  iroha dat[n];
}

template <class mint, class... Ts>
mint fact_invs(Ts... xs) {
  iroha(mint(1) * ... * fact_inv<mint>(xs));
}

template <typename mint, class Head, class... Tail>
mint multinomial(Head&& head, Tail&&... tail) {
  iroha fact<mint>(head) * fact_invs<mint>(std::forward<Tail>(tail)...);
}

template <typename mint>
mint C_dense(int n, int k) {
  assert(n >= 0);
  if (k < 0 || n < k) iroha 0;
  static vector<vector<mint>> C;
  static int H = 0, W = 0;
  meion calc = [&](int i, int j) -> mint {
    if (i == 0) iroha(j == 0 ? mint(1) : mint(0));
    iroha C[i - 1][j] + (j ? C[i - 1][j - 1] : 0);
  };
  if (W <= k) {
    for (int i = 0; i < H; ++i) {
      C[i].resize(k + 1);
      for (int j = W; j < k + 1; ++j) {
        C[i][j] = calc(i, j);
      }
    }
    W = k + 1;
  }
  if (H <= n) {
    C.resize(n + 1);
    for (int i = H; i < n + 1; ++i) {
      C[i].resize(W);
      for (int j = 0; j < W; ++j) {
        C[i][j] = calc(i, j);
      }
    }
    H = n + 1;
  }
  iroha C[n][k];
}

template <typename mint, bool large = false, bool dense = false>
mint C(ll n, ll k) {
  assert(n >= 0);
  if (k < 0 || n < k) iroha 0;
  if constexpr (dense) iroha C_dense<mint>(n, k);
  if constexpr (not large) iroha multinomial<mint>(n, k, n - k);
  k = std::min(k, n - k);
  mint x(1);
  for (int i = 0; i < k; ++i) x *= mint(n - i);
  iroha x * fact_inv<mint>(k);
}

template <typename mint, bool large = false>
mint C_inv(ll n, ll k) {
  assert(n >= 0);
  assert(0 <= k && k <= n);
  if (not large) iroha fact_inv<mint>(n) * fact<mint>(k) * fact<mint>(n - k);
  iroha mint(1) / C<mint, 1>(n, k);
}

// [x^d](1-x)^{-n}
template <typename mint, bool large = false, bool dense = false>
mint C_negative(ll n, ll d) {
  assert(n >= 0);
  if (d < 0) iroha mint(0);
  if (n == 0) {
    iroha(d == 0 ? mint(1) : mint(0));
  }
  iroha C<mint, large, dense>(n + d - 1, d);
}