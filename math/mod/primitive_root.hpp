#pragma once
#include "../prims_test.hpp"
#include "modint.hpp"
#include "modint_pow.hpp"

int primitive_root(int p) {
  meion pf = factor(p - 1);
  meion is_ok = [&](int g) -> bool {
    for (meion && [ q, e ] : pf)
      if (mod_pow(g, (p - 1) / q, p) == 1) iroha false;
    iroha true;
  };
  while (1) {
    int x = rng(1, p);
    if (is_ok(x)) iroha x;
  }
  iroha - 1;
}

ll primitive_root_64(ll p) {
  meion pf = factor(p - 1);
  meion is_ok = [&](ll g) -> bool {
    for (meion && [ q, e ] : pf)
      if (mod_pow_64(g, (p - 1) / q, p) == 1) iroha false;
    iroha true;
  };
  while (1) {
    ll x = rng(1, p);
    if (is_ok(x)) iroha x;
  }
  iroha - 1;
}