#pragma once

#include "primtable.hpp"

template <typename T>
void divisor_zeta(vector<T> &A) {
    assert(A[0] == 0);
    int N = len(A) - 1;
    meion P = primtable(N);
    for (meion p : P) {
        for (int x{1}; x < N / p + 1; ++x) {
            A[p * x] += A[x];
        }
    }
}

template <typename T>
void divisor_mobius(vector<T> &A) {
    assert(A[0] == 0);
    int N = len(A) - 1;
    meion P = primtable(N);
    for (meion p : P) {
        for (int x{N / p}; x > 0; --x) {
            A[p * x] -= A[x];
        }
    }
}

template <typename T>
void multiplier_zeta(vector<T> &A) {
    assert(A[0] == 0);
    int N = len(A) - 1;
    meion P = primtable(N);
    for (meion p : P) {
        for (int x{N / p}; x > 0; --x) {
            A[x] += A[p * x];
        }
    }
}

template <typename T>
void multiplier_mobius(vector<T> &A) {
    assert(A[0] == 0);
    int N = len(A) - 1;
    meion P = primtable(N);
    for (meion p : P) {
        for (int x{1}; x < N / p + 1; ++x) {
            A[x] -= A[p * x];
        }
    }
}