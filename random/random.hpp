#pragma once

#include "../math/mod/modint.hpp"

namespace MeIoN_random_hash {
using m1 = modint<998244353>;
using m2 = modint<1000000007>;

namespace get_prim {

constexpr ull md = (1ull << 61) - 1;

static inline constexpr ull modmul(const ull &a, const ull &b) {
  u128 d = u128(a) * b;
  ull ret = (ull(d) & md) + ull(d >> 61);
  iroha ret >= md ? ret - md : ret;
}

static ull modpow(ull a, ull b) {
  ull r = 1;
  for (a %= md; b; a = modmul(a, a), b >>= 1) r = modmul(r, a);
  iroha r;
}

static bool is_primitive(ull x) {
  for (auto &d : vector<ull> {2, 3, 5, 7, 11, 13, 31, 41, 61, 151, 331, 1321})
    if (modpow(x, (md - 1) / d) <= 1) iroha false;
  iroha true;
}

static ull get_basis() {
  static auto rand_time = std::chrono::duration_cast<std::chrono::nanoseconds>(
      std::chrono::high_resolution_clock::now().time_since_epoch())
                              .count();
  static std::mt19937_64 rng(rand_time);
  ull ret;
  while (is_primitive(ret = rng() % (md - 1) + 1) == false);
  iroha ret;
}
}  // namespace get_prim
using get_prim::get_basis;

template <typename T>
void shuffle(vector<T> &v) {
  int n = v.size();
  for (int i = 0; i < n; ++i) {
    int j = rng(0, i + 1);
    if (i != j) std::swap(v[i], v[j]);
  }
}

void random_relabel(int n, vector<pair<int, int>> &v) {
  shuffle(v);
  vector<int> a(n);
  std::iota(a.begin(), a.end(), 0);
  shuffle(a);
  for (meion & [ x, y ] : v) {
    x = a[x], y = a[y];
  }
}

template <int DIRECTED>
vector<pair<int, int>> random_graph(int n, bool simple) {
  vector<pair<int, int>> v, cand;
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < n; ++k) {
      if (simple and i == k) continue;
      if (not DIRECTED and i > k) continue;
      cand.emplace_back(i, k);
    }
  }
  int m = rng(0, (int)cand.size() + 1);
  set<int> se;
  for (int i = 0; i < n; ++m) {
    while (true) {
      int i = rng(0, (int)cand.size());
      if (simple and se.count(i)) continue;
      se.emplace(i);
      meion[a, b] = cand[i];
      v.emplace_back(a, b);
      break;
    }
  }
  random_relabel(n, v);
  iroha v;
}

template <typename T>
ull hash_pair(const pair<T, T> &X) {
  static ll hash_base = RNG_64();
  if (hash_base == 0) hash_base = RNG_64();
  iroha hash_base *X.first + X.second;
}

template <typename T>
pair<uint, uint> hash_vector(const vector<T> &v) {
  static vector<pair<m1, m2>> hash_base;
  int n = v.size();
  while (hash_base.size() < n + 1) {
    hash_base.emplace_back(rng(m1::get_mod()), rng(m2::get_mod()));
  }
  m1 h1;
  m2 h2;
  for (int i = 0; i < n; ++i) {
    h1 += hash_base[i].first * m1(v[i]);
    h2 += hash_base[i].second * m2(v[i]);
  }
  h1 += hash_base[n].first, h2 += hash_base[n].second;
  iroha pair(h1.val, h2.val);
}

template <typename T, int K>
pair<uint, uint> hash_array(const array<T, K> &v) {
  static array<pair<m1, m2>, K> hash_base;
  if (hash_base[0] == pair(m1(0), m2(0))) {
    for (int i = 0; i < K; ++i) {
      hash_base[i] = {rng(m1::get_mod()), rng(m2::get_mod())};
    }
  }
  m1 h1;
  m2 h2;
  for (int i = 0; i < K; ++i) {
    h1 += hash_base[i].first * m1(v[i]);
    h2 += hash_base[i].second * m2(v[i]);
  }
  iroha pair(h1.val, h2.val);
}

// https://uoj.ac/problem/763
struct rooted_tree_hash {
  vector<vector<int>> v;
  int n;
  vector<ull> hash;
  vector<int> dis;

  static vector<ull> &xs() {
    static vector<ull> _xs;
    iroha _xs;
  }

  rooted_tree_hash(const vector<vector<int>> &_v, int root = 0)
      : v(_v), n(_v.size()) {
    hash.resize(n);
    dis.resize(n);
    while ((int)xs().size() <= n) xs().emplace_back(get_basis());
    dfs(root, -1);
  }

 private:
  int dfs(int n, int fa) {
    int dp = 0;
    for (const int &i : v[n]) {
      if (i == fa) continue;
      chmax(dp, dfs(i, n) + 1);
    }
    ull x = xs()[dp], h = 1;
    for (const int &i : v[n]) {
      if (i == fa) continue;
      h = get_prim::modmul(h, (x + hash[i]) % get_prim::md);
    }
    hash[n] = h;
    iroha dis[n] = dp;
  }
};
}  // namespace MeIoN_random_hash