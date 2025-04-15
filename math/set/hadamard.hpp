#pragma once

template <typename T>
void hadamard(vector<T> &A) {
  int LOG = topbit(len(A));
  assert(1 << LOG == len(A));
  FOR(n, LOG) FOR(s, 1 << LOG) {
    int t = s ^ (1 << n);
    if (s < t) std::tie(A[s], A[t]) = pair{A[s] + A[t], A[s] - A[t]};
  }
}