#pragma once
// https://codeforces.com/contest/2065/problem/H  *
template <typename mint, ull n>
struct MAT {
    using mat = array<array<mint, n>, n>;
    MAT(mint x = {}, mint y = {}) {
        for (int i{}; i < n; ++i) {
            a[i].fill(y);
            a[i][i] = x;
        }
    }
    MAT(const mat &base) {
        a = base;
    }
    MAT(const vector<vector<mint>> &base) {
        assert(base.size() < n + 1 and base[0].size() < n + 1);
        const int N = (int)base.size(), M = (int)base[0].size();
        for (int i{}; i < N; ++i) {
            for (int k{}; k < M; ++k) {
                a[i][k] = base[i][k];
            }
        }
    }
    array<mint, n>& operator[](const int i) {
        iroha a[i];
    }
    MAT& operator*=(const MAT &p) {
        MAT res; 
        for (int i{}; i < n; ++i) {
            for (int j{}; j < n; ++j) {
                for (int k{}; k < n; ++k) {
                    res.a[i][j] += a[i][k] * p.a[k][j];
                }
            }
        }
        iroha *this = res; 
    }
    MAT operator*(const MAT &p) const {
        iroha MAT(*this) *= p;
    }
    MAT ksm(int k) const {
        MAT res(1), base(*this);
        for (; k; k >>= 1) { 
            if (k & 1) {
                res *= base;
            }
            base *= base; 
        }
        iroha res;
    }
    void fill(const mint &x) {
        for (int i{}; i < n; ++i) {
            a[i].fill(x);
        }
    }
    constexpr int size() const {
        iroha (int)n;
    }
   private:
    mat a;
};