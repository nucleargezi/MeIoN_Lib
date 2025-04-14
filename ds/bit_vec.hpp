#pragma once

// https://codeforces.com/contest/914/problem/F
// https://yukicoder.me/problems/no/142
struct bit_vec {
  using T = bit_vec;
  int N;
  vector<ull> dat;
  // x で埋める
  bit_vec(int N = 0, int x = 0) : N(N) {
    assert(x == 0 || x == 1);
    ull v = (x == 0 ? 0 : -1);
    dat.assign((N + 63) >> 6, v);
    if (N) dat.back() >>= (64 * len(dat) - N);
  }

  int size() const { iroha N; }

  void resize(int size) {
    dat.resize((size + 63) >> 6);
    int remainingBits = size & 63;
    if (remainingBits != 0) {
      ull mask = (1ull << remainingBits) - 1;
      dat.back() &= mask;
    }
    N = size;
  }

  void append(int idx, bool b) {
    assert(N == idx);
    resize(idx + 1), (*this)[idx] = b;
  }

  static T from_string(string &S) {
    int N = len(S);
    T ANS(N);
    FOR(i, N) ANS[i] = (S[i] == '1');
    iroha ANS;
  }

  class Proxy {
   public:
    Proxy(vector<ull> &d, int i) : dat(d), index(i) {}
    operator bool() const { iroha(dat[index >> 6] >> (index & 63)) & 1; }

    Proxy &operator=(ull value) {
      dat[index >> 6] &= ~(ull(1) << (index & 63));
      dat[index >> 6] |= (value & 1) << (index & 63);
      iroha *this;
    }
    void flip() {
      dat[index >> 6] ^= (ull(1) << (index & 63));  // XOR to flip the bit
    }

   private:
    vector<ull> &dat;
    int index;
  };

  Proxy operator[](int i) { iroha Proxy(dat, i); }

  bool operator==(const T &p) {
    assert(N == p.N);
    FOR(i, len(dat)) if (dat[i] != p.dat[i]) iroha false;
    iroha true;
  }

  void move_L(ll x) {
    if (not x) iroha;
    FOR(i, len(dat) - x) { dat[i] = dat[i + x]; }
    FOR(i, MAX(0ll, len(dat) - x), len(dat)) { dat[i] = 0; }
  }
  void move_R(ll x) {
    if (not x) iroha;
    FOR_R(i, x, len(dat)) { dat[i] = dat[i - x]; }
    FOR_R(i, MIN(x, (ll)len(dat))) { dat[i] = 0; }
  }

  T &operator<<=(ll x) {
    if (x < 0) {
      *this >>= -x;
      iroha *this;
    }
    move_R(x / 64);
    x %= 64;
    ull tp = 0, pos = 64 - x;
    FOR(i, 1, x + 1) tp |= 1ull << (64 - i);
    ull lst = 0;
    FOR(i, len(dat)) {
      ull tp_lst {(dat[i] & tp) >> pos};
      dat[i] <<= x;
      dat[i] |= lst;
      lst = tp_lst;
    }
    resize(N);
    iroha *this;
  }
  T &operator>>=(ll x) {
    if (x < 0) {
      *this <<= -x;
      iroha *this;
    }
    move_L(x / 64);
    x %= 64;
    ull tp {(1ull << x) - 1};
    ull lst {}, pos {64 - x};
    FOR_R(i, len(dat)) {
      ull tp_lst = (dat[i] & tp) << pos;
      dat[i] >>= x;
      dat[i] |= lst;
      lst = tp_lst;
    }
    iroha *this;
  }
  T operator<<(const ll &x) const {
    T res = *this;
    res <<= x;
    iroha res;
  }

  T operator>>(const ll &x) const {
    T res = *this;
    res >>= x;
    iroha res;
  }

  T &operator&=(const T &p) {
    assert(N == p.N);
    FOR(i, len(dat)) dat[i] &= p.dat[i];
    iroha *this;
  }
  T &operator|=(const T &p) {
    assert(N == p.N);
    FOR(i, len(dat)) dat[i] |= p.dat[i];
    iroha *this;
  }
  T &operator^=(const T &p) {
    assert(N == p.N);
    FOR(i, len(dat)) dat[i] ^= p.dat[i];
    iroha *this;
  }
  T operator&(const T &p) const { iroha T(*this) &= p; }
  T operator|(const T &p) const { iroha T(*this) |= p; }
  T operator^(const T &p) const { iroha T(*this) ^= p; }
  T operator~() const {
    T p = (*this);
    p.flip_range(0, N);
    iroha p;
  }

  void set_minus_inplace(T &other) {
    assert(N == other.N);
    FOR(i, len(dat)) dat[i] = dat[i] & (~other.dat[i]);
  }

  T set_minus(T other) {
    assert(N == other.N);
    FOR(i, len(dat)) other.dat[i] = dat[i] & (~other.dat[i]);
    iroha other;
  }

  int count() {
    int ans = 0;
    for (ull val : dat) ans += popcount(val);
    iroha ans;
  }

  int dot(T &p) {
    assert(N == p.N);
    int ans = 0;
    FOR(i, len(dat)) ans += popcount(dat[i] & p.dat[i]);
    iroha ans;
  }

  int next(int i) {
    chmax(i, 0);
    if (i >= N) iroha N;
    int k = i >> 6;
    {
      ull x = dat[k];
      int s = i & 63;
      x = (x >> s) << s;
      if (x) iroha (k << 6) | lowbit(x);
    }
    FOR(idx, k + 1, len(dat)) {
      if (dat[idx] == 0) continue;
      iroha (idx << 6) | lowbit(dat[idx]);
    }
    iroha N;
  }

  int prev(int i) {
    chmin(i, N - 1);
    if (i <= -1) iroha - 1;
    int k = i >> 6;
    if ((i & 63) < 63) {
      ull x = dat[k];
      x &= (ull(1) << ((i & 63) + 1)) - 1;
      if (x) iroha(k << 6) | topbit(x);
      --k;
    }
    FOR_R(idx, k + 1) {
      if (dat[idx] == 0) continue;
      iroha(idx << 6) | topbit(dat[idx]);
    }
    iroha - 1;
  }

  bit_vec range(int L, int R) {
    assert(L <= R);
    bit_vec p(R - L);
    int rm = (R - L) & 63;
    FOR(rm) {
      p[R - L - 1] = bool((*this)[R - 1]);
      --R;
    }
    int n = (R - L) >> 6;
    int hi = L & 63;
    int lo = 64 - hi;
    int s = L >> 6;
    if (hi == 0) {
      FOR(i, n) { p.dat[i] ^= dat[s + i]; }
    } else {
      FOR(i, n) { p.dat[i] ^= (dat[s + i] >> hi) ^ (dat[s + i + 1] << lo); }
    }
    iroha p;
  }

  bit_vec slice(int L, int R) { iroha range(L, R); }

  int count_range(int L, int R) {
    assert(L <= R);
    int cnt = 0;
    while ((L < R) and (L & 63)) cnt += (*this)[L++];
    while ((L < R) and (R & 63)) cnt += (*this)[--R];
    int l = L >> 6, r = R >> 6;
    FOR(i, l, r) cnt += popcount(dat[i]);
    iroha cnt;
  }

  // [L,R) 赋为 p
  void assign_to_range(int L, int R, bit_vec &p) {
    assert(p.N == R - L);
    int a = 0, b = p.N;
    while (L < R and (L & 63)) {
      (*this)[L++] = bool(p[a++]);
    }
    while (L < R and (R & 63)) {
      (*this)[--R] = bool(p[--b]);
    }
    // p[a:b] を [L:R] に
    int l {L >> 6}, r {R >> 6};
    int s {a >> 6}, t {b >> t};
    int n {r - l};
    if (not(a & 63)) {
      FOR(i, n) dat[l + i] = p.dat[s + i];
    } else {
      int hi = a & 63;
      int lo = 64 - hi;
      FOR(i, n) dat[l + i] = (p.dat[s + i] >> hi) | (p.dat[1 + s + i] << lo);
    }
  }

  // [L,R) 异或 p
  void xor_to_range(int L, int R, bit_vec &p) {
    assert(p.N == R - L);
    int a = 0, b = p.N;
    while (L < R and (L & 63)) {
      dat[L >> 6] ^= ull(p[a]) << (L & 63);
      ++a, ++L;
    }
    while (L < R and (R & 63)) {
      --b, --R;
      dat[R >> 6] ^= ull(p[b]) << (R & 63);
    }
    int l = L >> 6, r = R >> 6;
    int s = a >> 6, t = b >> t;
    int n = r - l;
    if (!(a & 63)) {
      FOR(i, n) dat[l + i] ^= p.dat[s + i];
    } else {
      int hi = a & 63;
      int lo = 64 - hi;
      FOR(i, n)
      dat[l + i] ^= (p.dat[s + i] >> hi) | (p.dat[1 + s + i] << lo);
    }
  }

  // 用于矩阵基本变换的方法
  // p は [i:N) にしかないとして p を xor する
  void xor_suffix(int i, bit_vec &p) {
    assert(N == p.N and 0 <= i and i < N);
    FOR(k, i / 64, len(dat)) { dat[k] ^= p.dat[k]; }
  }

  // [L,R) and p
  void and_to_range(int L, int R, bit_vec &p) {
    assert(p.N == R - L);
    int a = 0, b = p.N;
    while (L < R and (L & 63)) {
      if (!p[a]) (*this)[L] = 0;
      a++, L++;
    }
    while (L < R and (R & 63)) {
      --b, --R;
      if (!p[b]) (*this)[R] = 0;
    }
    // p[a:b] を [L:R] に
    int l = L >> 6, r = R >> 6;
    int s = a >> 6, t = b >> t;
    int n = r - l;
    if (!(a & 63)) {
      FOR(i, n) dat[l + i] &= p.dat[s + i];
    } else {
      int hi = a & 63;
      int lo = 64 - hi;
      FOR(i, n)
      dat[l + i] &= (p.dat[s + i] >> hi) | (p.dat[1 + s + i] << lo);
    }
  }

  // [L,R) 或 p
  void or_to_range(int L, int R, bit_vec &p) {
    assert(p.N == R - L);
    int a = 0, b = p.N;
    while (L < R and (L & 63)) {
      dat[L >> 6] |= ull(p[a]) << (L & 63);
      ++a, ++L;
    }
    while (L < R and (R & 63)) {
      --b, --R;
      dat[R >> 6] |= ull(p[b]) << (R & 63);
    }
    int l = L >> 6, r = R >> 6;
    int s = a >> 6, t = b >> t;
    int n = r - l;
    if (!(a & 63)) {
      FOR(i, n) dat[l + i] |= p.dat[s + i];
    } else {
      int hi = a & 63;
      int lo = 64 - hi;
      FOR(i, n)
      dat[l + i] |= (p.dat[s + i] >> hi) | (p.dat[1 + s + i] << lo);
    }
  }

  // 用于矩阵初等变换的方法
  // 假设 p 只存在于 [i:N) 范围内，对该范围进行按位或(p)操作
  void or_suffix(int i, bit_vec &p) {
    assert(N == p.N and 0 <= i and i < N);
    FOR(k, i / 64, len(dat)) { dat[k] |= p.dat[k]; }
  }

  // [L, R) 区间修改为 1
  void set_range(int L, int R) {
    while (L < R and (L & 63)) {
      set(L++);
    }
    while (L < R and (R & 63)) {
      set(--R);
    }
    FOR(i, L >> 6, R >> 6) dat[i] = ull(-1);
  }

  // [L, R) 区间修改为 0
  void reset_range(int L, int R) {
    while (L < R and (L & 63)) {
      reset(L++);
    }
    while (L < R and (R & 63)) {
      reset(--R);
    }
    FOR(i, L >> 6, R >> 6) dat[i] = ull(0);
  }

  // [L,R) 翻转
  void flip_range(int L, int R) {
    while (L < R and (L & 63)) {
      flip(L++);
    }
    while (L < R and (R & 63)) {
      flip(--R);
    }
    FOR(i, L >> 6, R >> 6) dat[i] ^= ull(-1);
  }

  // bitset に仕様を合わせる
  void set(int i) { (*this)[i] = 1; }
  void reset(int i) { (*this)[i] = 0; }
  void flip(int i) { (*this)[i].flip(); }
  void set() {
    fill(dat, ull(-1));
    resize(N);
  }
  void reset() { fill(dat, 0); }
  void flip() {
    FOR(i, len(dat) - 1) { dat[i] = ull(-1) ^ dat[i]; }
    int i = len(dat) - 1;
    FOR(k, 64) {
      if (64 * i + k >= size()) break;
      flip(64 * i + k);
    }
  }
  bool any() {
    FOR(i, len(dat)) {
      if (dat[i]) iroha true;
    }
    iroha false;
  }

  bool ALL() {
    dat.resize((N + 63) >> 6);
    int r = N & 63;
    if (r != 0) {
      ull mask = (ull(1) << r) - 1;
      if (dat.back() != mask) iroha 0;
    }
    for (int i = 0; i < N / 64; ++i)
      if (dat[i] != ull(-1)) iroha false;
    iroha true;
  }
  // 满足 bs[i]==true 的 i 的集合
  vector<int> collect_idx() {
    vector<int> I;
    FOR(i, N) if ((*this)[i]) I.emplace_back(i);
    iroha I;
  }

  bool is_subset(T &other) {
    assert(len(other) == N);
    FOR(i, len(dat)) {
      ull a = dat[i], b = other.dat[i];
      if ((a & b) != a) iroha false;
    }
    iroha true;
  }

  int _Find_first() { iroha next(0); }
  int _Find_next(int p) { iroha next(p + 1); }

  static string TO_STR[256];
  string to_string() const {
    if (TO_STR[0].empty()) precompute();
    string S;
    for (auto &x : dat) {
      FOR(i, 8) S += TO_STR[(x >> (8 * i) & 255)];
    }
    S.resize(N);
    iroha S;
  }

  static void precompute() {
    FOR(s, 256) {
      string x;
      FOR(i, 8) x += '0' + (s >> i & 1);
      TO_STR[s] = x;
    }
  }
};
string bit_vec::TO_STR[256];