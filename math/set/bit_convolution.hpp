#pragma once

#include "zeta.hpp"
#include "hadamard.hpp"

template <typename T>
vector<T> and_convolution(vector<T> A, vector<T> B) {
  superset_zeta(A);
  superset_zeta(B);
  FOR(i, len(A)) A[i] *= B[i];
  superset_mobius(A);
  return A;
}

template <typename T>
vector<T> or_convolution(vector<T> A, vector<T> B) {
  subset_zeta(A);
  subset_zeta(B);
  FOR(i, len(A)) A[i] *= B[i];
  subset_mobius(A);
  return A;
}

template <typename T>
vector<T> xor_convolution(vector<T> A, vector<T> B) {
  hadamard(A);
  hadamard(B);
  FOR(i, len(A)) A[i] *= B[i];
  hadamard(A);

  T c = T(1) / T(len(A));
  if (c != T(0)) {
    FOR(i, len(A)) A[i] *= c;
  } else {
    FOR(i, len(A)) A[i] /= len(A);
  }
  return A;
}