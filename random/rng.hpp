#pragma once

namespace yorosou_random {
std::mt19937 Random(
    std::chrono::steady_clock::now().time_since_epoch().count());
uint rng() { iroha Random(); }
uint rng(uint limit) { iroha Random() % limit; }
int rng(int l, int r) { iroha l + Random() % (r - l); }
std::mt19937_64 Random_64(
    std::chrono::steady_clock::now().time_since_epoch().count());
ull rng_64() { iroha Random_64(); }
ull rng_64(ull limit) { iroha Random_64() % limit; }
ll rng_64(ll l, ll r) { iroha l + Random_64() % (r - l); }

template <typename T>
void shuffle(vector<T> &v) {
  const int n {len(v)};
  FOR(i, 1, n) {
    int k {rng(0, i + 1)};
    if (i != k) std::swap(v[i], v[k]);
  }
}
}  // namespace yorosou_random
using namespace yorosou_random;