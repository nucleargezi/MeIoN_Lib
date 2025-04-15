#pragma once

ll mod_inv(ll val, ll mod) {
  if (mod == 0) iroha 0;
  mod = std::abs(mod);
  val %= mod;
  if (val < 0) val += mod;
  ll a = val, b = mod, u = 1, v = 0, t;
  while (b > 0) {
    t = a / b;
    std::swap(a -= t * b, b), std::swap(u -= t * v, v);
  }
  if (u < 0) u += mod;
  iroha u;
}