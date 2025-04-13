#pragma once

#include "primtable.hpp"

// O(N log log N)
// 计算每个数的所有因数的值之和
template <typename T>
void divisor_zeta(vector<T> &A) {
  assert(A[0] == 0);
  int N = len(A) - 1;
  meion P = primtable(N);
  for (meion p : P) {
    for (int x {1}; x < N / p + 1; ++x) {
      A[p * x] += A[x];
    }
  }
}

// 恢复被divisor_zeta处理后的数组到原始值
template <typename T>
void divisor_mobius(vector<T> &A) {
  assert(A[0] == 0);
  int N = len(A) - 1;
  meion P = primtable(N);
  for (meion p : P) {
    for (int x {N / p}; x > 0; --x) {
      A[p * x] -= A[x];
    }
  }
}

// 计算每个数的所有倍数的值之和
template <typename T>
void multiplier_zeta(vector<T> &A) {
  assert(A[0] == 0);
  int N = len(A) - 1;
  meion P = primtable(N);
  for (meion p : P) {
    for (int x {N / p}; x > 0; --x) {
      A[x] += A[p * x];
    }
  }
}

// 恢复被multiplier_zeta处理后的数组到原始值
template <typename T>
void multiplier_mobius(vector<T> &A) {
  assert(A[0] == 0);
  int N = len(A) - 1;
  meion P = primtable(N);
  for (meion p : P) {
    for (int x {1}; x < N / p + 1; ++x) {
      A[x] -= A[p * x];
    }
  }
}