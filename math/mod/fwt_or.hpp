#pragma once

template <typename mint>
void fwt_or(vector<mint> &a) {
    const int n = (int)a.size();
    for (int msk = 1; msk < n; msk <<= 1) {
        for (int i{}; i < n; i += msk << 1) {
            for (int k{}; k < msk; k++) {
                a[msk + i + k] += a[i + k];
            }
        }
    }
}
template <typename mint>
void fwt_ior(vector<mint> &a) {
    const int n = (int)a.size();
    for (int bit{1}; bit < n; bit <<= 1) {
        for (int i{}; i < n; i += bit << 1) {
            for (int k{}; k < bit; k++) {
                a[bit + i + k] -= a[i + k];
            }
        }
    }
}