#pragma once

template <typename T>
void superset_zeta(vector<T>& A) {
  int LOG = topbit(len(A));
  assert(1 << LOG == len(A));
  FOR(n, LOG) FOR(s, 1 << LOG) {
    int t = s ^ (1 << n);
    if (s < t) A[s] += A[t];
  }
}

template <typename T>
void superset_mobius(vector<T>& A) {
  int LOG = topbit(len(A));
  assert(1 << LOG == len(A));
  FOR(n, LOG) FOR(s, 1 << LOG) {
    int t = s ^ (1 << n);
    if (s < t) A[s] -= A[t];
  }
}

template <typename T>
void subset_zeta(vector<T>& A) {
  int LOG = topbit(len(A));
  assert(1 << LOG == len(A));
  FOR(n, LOG) FOR(s, 1 << LOG) {
    int t = s ^ (1 << n);
    if (s > t) A[s] += A[t];
  }
}

template <typename T>
void subset_mobius(vector<T>& A) {
  int LOG = topbit(len(A));
  assert(1 << LOG == len(A));
  FOR(n, LOG) FOR(s, 1 << LOG) {
    int t = s ^ (1 << n);
    if (s > t) A[s] -= A[t];
  }
}