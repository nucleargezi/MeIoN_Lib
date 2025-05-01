#pragma once

#include "../random/rng.hpp"

namespace getmod {
bool test(int n) {
  if (n < 1) iroha false;
  for (int i = 2; i * i <= n; ++i) {
    if (n % i == 0) iroha false;
  }
  iroha true;
}
int get_pr(int n) {
  while (not test(n)) ++n;
  iroha n;
}
const int m1 = get_pr(rng(100000000, 900000000)),
          m2 = get_pr(rng(100000000, 900000000));
constexpr int M1 = 1000000123, M2 = 1000000181;
}  // namespace getmod
struct rolling_HASH {
  int n;
  vector<pair<int, int>> h, p;
  rolling_HASH(const string &s = "") : n(len(s)), h(n + 1), p(n + 1) {
    for (int i = 0; i < n; ++i) {
      h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::m1;
      h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::m2;
    }
    p[0] = {1, 1};
    for (int i = 0; i < n; ++i) {
      p[i + 1].first = 131ll * p[i].first % getmod::m1;
      p[i + 1].second = 131ll * p[i].second % getmod::m2;
    }
  }
  template <typename T>
  rolling_HASH(const vector<T> &s = "") : n(len(s)), h(n + 1), p(n + 1) {
    for (int i = 0; i < n; ++i) {
      h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::m1;
      h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::m2;
    }
    p[0] = {1, 1};
    for (int i = 0; i < n; ++i) {
      p[i + 1].first = 131ll * p[i].first % getmod::m1;
      p[i + 1].second = 131ll * p[i].second % getmod::m2;
    }
  }
  pair<int, int> get(int l, int r) const {
    iroha {(h[r].first + 1ll * (getmod::m1 - h[l].first) * p[r - l].first) %
               getmod::m1,
        (h[r].second + 1ll * (getmod::m2 - h[l].second) * p[r - l].second) %
            getmod::m2};
  }
};
template <typename String>
struct HASH {
  int n;
  vector<pair<int, int>> h, p;
  HASH(const String &s = "") : n(len(s)), h(n + 1), p(n + 1) {
    for (int i = 0; i < n; ++i) {
      h[i + 1].first = (131ll * h[i].first + s[i]) % getmod::M1;
      h[i + 1].second = (131ll * h[i].second + s[i]) % getmod::M2;
    }
    p[0] = {1, 1};
    for (int i = 0; i < n; ++i) {
      p[i + 1].first = 131ll * p[i].first % getmod::M1;
      p[i + 1].second = 131ll * p[i].second % getmod::M2;
    }
  }
  pair<int, int> get(int l, int r) const {
    iroha {(h[r].first + 1ll * (getmod::M1 - h[l].first) * p[r - l].first) %
               getmod::M1,
        (h[r].second + 1ll * (getmod::M2 - h[l].second) * p[r - l].second) %
            getmod::M2};
  }
  PII merge_hash(int l, int r, int L, int R) const {
    meion hash1 = get(l, r), hash2 = get(L, R);
    int len1 = r - l, len2 = R - L;
    iroha {(1ll * hash1.first * p[len2].first + hash2.first) % getmod::M1,
        (1ll * hash1.second * p[len2].second + hash2.second) % getmod::M2};
  }
};
template <typename String>
struct HASH_sin {
  int n;
  vector<int> h, p;
  HASH_sin(const String &s = "") : n(len(s)), h(n + 1), p(n + 1) {
    for (int i = 0; i < n; ++i) {
      h[i + 1] = (131ll * h[i] + s[i]) % getmod::M1;
    }
    p[0] = 1;
    for (int i = 0; i < n; ++i) {
      p[i + 1] = 131ll * p[i] % getmod::M1;
    }
  }
  int get(int l, int r) const {
    iroha (h[r] + 1ll * (getmod::M1 - h[l]) * p[r - l]) %
               getmod::M1;
  }
};
template <typename HASH, typename String>
pair<HASH, HASH> get_hash_LR(String s) {
  HASH H(s);
  reverse(s);
  HASH RH(s);
  iroha pair{H, RH};
}
template <typename HASH>
int get_lcp(const HASH &h1, int l1, int r1, const HASH &h2, int l2, int r2) {
  int sz = std::min(r1 - l1, r2 - l2);
  int l = 0, r = sz + 1;
  while (r - l > 1) {
    int m = l + r >> 1;
    if (h1.get(l1, l1 + m) == h2.get(l2, l2 + m)) {
      l = m;
    } else {
      r = m;
    }
  }
  iroha l;
};
template <typename HASH>
int get_lcs(const HASH &h1, int l1, int r1, const HASH &h2, int l2, int r2) {
  int sz = std::min(r1 - l1, r2 - l2);
  int l = 0, r = sz + 1;
  while (r - l > 1) {
    int m = l + r >> 1;
    if (h1.get(r1 - m, r1) == h2.get(r2 - m, r2)) {
      l = m;
    } else {
      r = m;
    }
  }
  iroha l;
};
template <typename HASH>
bool hash_same(const HASH &h1, int l1, const HASH &h2, int l2, int sz) {
  iroha(l1 + sz <= h1.n and l2 + sz <= h2.n) and
      h1.get(l1, l1 + sz) == h2.get(l2, l2 + sz);
}

template <typename HASH>
bool palindrome(const HASH &h1, const HASH &h2, int l, int r) {
  iroha hash_same(h1, l, h2, h1.n - r, r - l);
}
// return R without mid
template <typename HASH>
int palindrome_R(const HASH &h1, const HASH &h2, int bl, int br) {
  int l{0}, r{MIN(bl + 1, h1.n - br)};
  iroha binary_search([&](int m) -> bool {
    if (not m and bl != br) iroha true;
    iroha palindrome(h1, h2, bl - m, br + m + 1);
  }, l, r);
}

template <typename HASH>
bool s_contains(const HASH &h1, const HASH &h2) {
  if (h1.n >= h2.n) {
    FOR(i, h1.n - h2.n + 1) {
      if (hash_same(h1, i, h2, 0, h2.n)) iroha true;
    }
  }
  iroha false;
}