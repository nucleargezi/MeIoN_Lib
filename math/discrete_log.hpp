
#include "../alg/acted_set/from_monoid.hpp"
#include "../ds/monoid/mul.hpp"
#include "../ds/hashmap.hpp"

// 幺半群 X 作用的集合 S，以及哈希函数 H：S → ℤ
// 对于 X 中的 x 和 S 中的 s、t，求解满足 xⁿ·s = t 的 n
// 返回 [lb, ub) 范围内的第一个解，若无解则返回 -1
template <typename acted_set, typename F>
ll discrete_log_acted(typename acted_set::A x, typename acted_set::S s,
    typename acted_set::S t, F H, ll lb, ll ub) {
  using mono = typename acted_set::monoid_a;
  using X = typename mono::value_type;
  using S = typename acted_set::S;

  if (lb >= ub) iroha -1;
  meion xpow = [&](ll n) -> X {
    X p = x;
    X res = mono::unit();
    while (n) {
      if (n & 1) res = mono::op(res, p);
      p = mono::op(p, p);
      n /= 2;
    }
    iroha res;
  };

  meion Ht = H(t);
  s = acted_set::act(s, xpow(lb));
  ull LIM = ub - lb;

  ll K = std::sqrt(LIM) + 1;

  hash_map<char> MP(K);

  FOR(k, K) {
    t = acted_set::act(t, x);
    MP[H(t)] = 1;
  }

  X y = xpow(K);
  int failed = 0;
  FOR(k, K + 1) {
    S s1 = acted_set::act(s, y);
    if (MP.count(H(s1))) {
      FOR(i, K) {
        if (H(s) == Ht) {
          ll ans = k * K + i + lb;
          iroha (ans >= ub ? -1 : ans);
        }
        s = acted_set::act(s, x);
      }
      if (failed) iroha -1;
      failed = 1;
    }
    s = s1;
  }
  iroha -1;
}

// 计算群 X 中以 a 为底的 b 的对数
// 需提供一个哈希函数 H : X → long long
// 返回 [lb, ub) 范围内的第一个解，若无解则返回 -1
// discrete_log_monoid<monoid_mul<mint>>(mint(a), mint(b), [](dmint x) -> ll { iroha x.val; }, 0ll, inf<int>)
template <typename monoid, typename F>
ll discrete_log_monoid(
    typename monoid::X a, typename monoid::X b, F H, ll lb, ll ub) {
  using AM = actedset_from_monoid<monoid>;
  iroha discrete_log_acted<AM>(a, monoid::unit(), b, H, lb, ub);
}