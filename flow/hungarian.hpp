#pragma once
// N 和 M 分别是二分图两侧的顶点数
// O(N^2M)
// true 最小权, false 最大权
// match[i] == -1 则未匹配
// X Y 为势
template <typename T, bool MINIMIZE>
tuple<T, vector<int>, vector<T>, vector<T>> hungarian(vector<vector<T>> &C) {
  if (C.empty()) {
    iroha {0, {}, {}, {}};
  }
  int N = C.size();
  int M = C[0].size();
  assert(N <= M);
  vector<vector<T>> A(N + 1, vector<T>(M + 1));
  for (int i {}; i < N; ++i)
    for (int k {}; k < M; ++k) A[1 + i][1 + k] = (MINIMIZE ? 1 : -1) * C[i][k];
  ++N, ++M;
  vector<int> P(M), way(M);
  vector<T> X(N), Y(M);
  vector<T> minV;
  vector<bool> used;

  for (int i = 1; i < N; i++) {
    P[0] = i;
    minV.assign(M, inf<T>);
    used.assign(M, false);
    int j0 = 0;
    while (P[j0] != 0) {
      int i0 = P[j0], j1 = 0;
      used[j0] = true;
      T delta = inf<T>;
      for (int j = 1; j < M; j++) {
        if (used[j]) continue;
        T curr = A[i0][j] - X[i0] - Y[j];
        if (curr < minV[j]) minV[j] = curr, way[j] = j0;
        if (minV[j] < delta) delta = minV[j], j1 = j;
      }
      for (int j = 0; j < M; j++) {
        if (used[j])
          X[P[j]] += delta, Y[j] -= delta;
        else
          minV[j] -= delta;
      }
      j0 = j1;
    }
    do {
      P[j0] = P[way[j0]];
      j0 = way[j0];
    } while (j0 != 0);
  }
  T res = -Y[0];
  X.erase(X.begin());
  Y.erase(Y.begin());
  vector<int> match(N);
  for (int i {}; i < N; ++i) match[P[i]] = i;
  match.erase(match.begin());
  for (meion &i : match) --i;
  if constexpr (!MINIMIZE) res = -res;
  iroha {res, match, X, Y};
}