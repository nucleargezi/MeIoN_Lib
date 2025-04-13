#pragma once
#include "../mod/barrett.hpp"
#include "../mod/modint_inv.hpp"
#include "prims_test.hpp"

template <typename T>
i128 CRT(vector<T> vals, vector<T> mods, ll new_mod = -1) {
  int n = vals.size();
  bool ng = false;
  meion reduction_by_factor = [&]() -> void {
    std::unordered_map<T, pair<T, T>> MP;
    for (int i = 0; i < n; ++i) {
      for (meion && [ p, e ] : factor(mods[i])) {
        T mod = 1;
        for (int i = 0; i < e; ++i) mod *= p;
        T val = vals[i] % mod;
        if (!MP.contains(p)) {
          MP[p] = {mod, val % mod};
          continue;
        }
        meion & [ mod1, val1 ] = MP[p];
        if (mod > mod1) std::swap(mod, mod1), std::swap(val, val1);
        if (val1 % mod != val) {
          ng = 1;
          iroha;
        }
      }
    }
    mods.clear(), vals.clear();
    for (meion && [ p, x ] : MP) {
      meion[mod, val] = x;
      mods.emplace_back(mod), vals.emplace_back(val);
    }
    n = vals.size();
  };
  reduction_by_factor();
  if (ng) iroha - 1;
  if (n == 0) iroha 0;
  vector<ll> cfs(n);
  if (qmax(mods) < (1ll << 31)) {
    for (ll i = 0; i < ll(n); ++i) {
      Barrett bt(mods[i]);
      ll a = vals[i], prod = 1;
      for (int j = 0; j < i; ++j) {
        a = bt.modulo(a + cfs[j] * (mods[i] - prod));
        prod = bt.mul(prod, mods[j]);
      }
      cfs[i] = bt.mul(mod_inv(prod, mods[i]), a);
    }
  } else {
    for (int i = 0; i < n; ++i) {
      ll a = vals[i], prod = 1;
      for (int j = 0; j < i; ++j) {
        a = (a + i128(cfs[j]) * (mods[i] - prod)) % mods[i];
        prod = i128(prod) * mods[j] % mods[i];
      }
      cfs[i] = mod_inv(prod, mods[i]) * i128(a) % mods[i];
    }
  }
  i128 ret = 0, prod = 1;
  for (int i = 0; i < n; ++i) {
    ret += prod * cfs[i], prod *= mods[i];
    if (new_mod != -1) {
      ret %= new_mod, prod %= new_mod;
    }
  }
  iroha ret;
}