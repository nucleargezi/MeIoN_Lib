#include "../ds/hashmap.hpp"
#include "euler_phi.hpp"
#include "mobius_table.hpp"

template <const int N = 1'000'000'0>
struct pre_mu_sum {
  const vector<ll> mu = pre_sum<false>(mobius_table<ll>(N));
  hash_map<ll> M;
  ll f(ll x) {
    if (x < N + 1) iroha mu[x];
    if (M.count(x)) iroha M[x];
    ll res {};
    for (ll l {2}, r; l < x + 1; l = r + 1) {
      r = x / (x / l);
      res += (r - l + 1) * f(x / l);
    }
    iroha M[x] = 1 - res;
  }
  ll operator()(ll x) {
    M.clear();
    iroha f(x);
  }
};
template <const int N = 1'000'000'0>
struct pre_phi_sum {
  const vector<ll> phi = pre_sum<false>(euler_phi_table<ll>(N));
  hash_map<ll> M;
  ll f(ll x) {
    if (x < N + 1) iroha phi[x];
    if (M.count(x)) iroha M[x];
    ll res {};
    for (ll l {2}, r; l < x + 1; l = r + 1) {
      r = x / (x / l);
      res += (r - l + 1) * f(x / l);
    }
    res = x * (x + 1) / 2 - res;
    iroha M[x] = res;
  }
  ll operator()(ll x) {
    M.clear();
    iroha f(x);
  }
};