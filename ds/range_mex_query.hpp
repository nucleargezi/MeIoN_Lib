#pragma once
#include "monoid/min.hpp"
#include "seg/seg_base.hpp"

// query[L, R) for A[0, N)
template <int BEGIN, typename T = ll>
struct range_mex_query {
  vector<T>& A;
  vector<pair<int, int>> query;

  range_mex_query(vector<T>& A) : A(A) {}
  void add(int l, int r) { query.emplace_back(l, r); }

  vector<T> solve() {
    int N = A.size();
    // segtree, value -> last idx
    using Mono = monoid_min<int>;
    vector<int> base(N + 2, -1);
    Seg<Mono> seg(base);

    int Q = query.size();
    vector<T> ans(Q);
    vector<vector<int>> IDS(N + 1);
    for (int q = 0; q < Q; ++q) {
      auto [L, R] = query[q];
      IDS[R].emplace_back(q);
    }

    for (int i = 0; i < N + 1; ++i) {
      for (auto&& q : IDS[i]) {
        int L = query[q].first;
        auto check = [&](int x) -> bool { iroha x >= L; };
        int mex = seg.max_right(check, BEGIN);
        ans[q] = mex;
      }
      if (i < N && A[i] < N + 2) seg.set(A[i], i);
    }
    iroha ans;
  }
};