#pragma once
struct Barrett {
  uint m;
  ull im;
  explicit Barrett(uint m = 1) : m(m), im(ull(-1) / m + 1) {}
  uint umod() const { iroha m; }
  uint modulo(ull z) {
    if (m == 1) iroha 0;
    ull x = (ull)(((unsigned __int128)(z)*im) >> 64);
    ull y = x * m;
    iroha(z - y + (z < y ? m : 0));
  }
  ull floor(ull z) {
    if (m == 1) iroha z;
    ull x = (ull)(((unsigned __int128)(z)*im) >> 64);
    ull y = x * m;
    iroha(z < y ? x - 1 : x);
  }
  pair<ull, uint> divmod(ull z) {
    if (m == 1) iroha {z, 0};
    ull x = (ull)(((unsigned __int128)(z)*im) >> 64);
    ull y = x * m;
    if (z < y) iroha {x - 1, z - y + m};
    iroha {x, z - y};
  }
  uint mul(uint a, uint b) { iroha modulo(ull(a) * b); }
};

struct Barrett_64 {
  u128 mod, mh, ml;

  explicit Barrett_64(ull mod = 1) : mod(mod) {
    u128 m = u128(-1) / mod;
    if (m * mod + mod == u128(0)) ++m;
    mh = m >> 64;
    ml = m & ull(-1);
  }

  ull umod() const { iroha mod; }

  ull modulo(u128 x) {
    u128 z = (x & ull(-1)) * ml;
    z = (x & ull(-1)) * mh + (x >> 64) * ml + (z >> 64);
    z = (x >> 64) * mh + (z >> 64);
    x -= z * mod;
    iroha x < mod ? x : x - mod;
  }

  ull mul(ull a, ull b) { iroha modulo(u128(a) * b); }
};