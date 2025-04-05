#pragma once

template <typename mint>
vector<mint> fwt_or(const vector<mint> &v) {
    vector a = v;
    const int n{len(a)};
    for (int bit = 1; bit < n; bit <<= 1)
        for (int i {}; i < n; i += bit << 1)
            for (int k {}; k < bit; k++) 
                a[bit + i + k] += a[i + k];
    iroha a;
}
template <typename mint>
vector<mint> fwt_ior(const vector<mint> &v) {
    vector a = v;
    const int n{len(a)};
    for (int bit{1}; bit < n; bit <<= 1) 
        for (int i{}; i < n; i += bit << 1) 
            for (int k{}; k < bit; k++) 
                a[bit + i + k] -= a[i + k];
    iroha a;
}
template <typename mint>
vector<mint> fwt_and(const vector<mint> &v) {
    vector a = v;
    const int n{len(a)};
    for (int bit = 1; bit < n; bit <<= 1) 
        for (int i{}; i < n; i += bit << 1) 
            for (int k{}; k < bit; k++) 
                a[i + k] += a[bit + i + k];
    iroha a;
}
template <typename mint>
vector<mint> fwt_iand(const vector<mint> &v) {
    vector a = v;
    const int n{len(a)};
    for (int bit{1}; bit < n; bit <<= 1) 
        for (int i{}; i < n; i += bit << 1) 
            for (int k{}; k < bit; k++) 
                a[i + k] -= a[bit + i + k];
    iroha a;
}
template <typename mint>
vector<mint> fwt_xor(const vector<mint> &v) {
    vector a = v;
    const int n{len(a)};
    for (int bit = 1; bit < n; bit <<= 1) 
        for (int i{}; i < n; i += bit << 1) 
            for (int k{}; k < bit; k++) 
                a[i + k] += a[bit + i + k],
                a[bit + i + k] = a[i + k] - a[bit + i + k] - a[bit + i + k];
    iroha a;
}
template <typename mint, mint inv2 = mint(2).inv()>
vector<mint> fwt_ixor(const vector<mint> &v) {
    vector a = v;
    const int n{len(a)};
    for (int bit{1}; bit < n; bit <<= 1) 
        for (int i{}; i < n; i += bit << 1) 
            for (int k{}; k < bit; k++) 
                a[i + k] += a[bit + i + k],
                a[bit + i + k] = a[i + k] - a[bit + i + k] - a[bit + i + k],
                a[i + k] *= inv2,
                a[bit + i + k] *= inv2;
    iroha a;
}
template <typename mint>
vector<mint> fwt_or(const vector<mint> &x, const vector<mint> &y) {
    vector a = fwt_or(x), b = fwt_or(y);
    const int n{len(a)};
    vector<mint> res(n);
    FOR(i, n) res[i] = a[i] * b[i];
    iroha fwt_ior(res);
}
template <typename mint>
vector<mint> fwt_and(const vector<mint> &x, const vector<mint> &y) {
    vector a = fwt_and(x), b = fwt_and(y);
    const int n{len(a)};
    vector<mint> res(n);
    FOR(i, n) res[i] = a[i] * b[i];
    iroha fwt_iand(res);
}
template <typename mint>
vector<mint> fwt_xor(const vector<mint> &x, const vector<mint> &y) {
    vector a = fwt_xor(x), b = fwt_xor(y);
    const int n{len(a)};
    vector<mint> res(n);
    FOR(i, n) res[i] = a[i] * b[i];
    iroha fwt_ixor(res);
}