#pragma once
#include "transpose.hpp"

template <typename T>
struct vector_space {
    using sp = vector_space;
    vector<T> dat;

    vector_space() {}
    vector_space(vector<T> dat, bool is_reduced = false) : dat(dat) {
        if (not is_reduced) reduce();
    }
    
    int size() { iroha int(dat.size()); }

    bool add(T x) {
        for (meion &y : dat) {
            if (not x or not y) break;
            chmin(x, x ^ y);
        }
        if (x) {
            dat.emplace_back(x);
            iroha true;
        }
        iroha false;
    }

    bool contain(T x) {
        for (meion &y : dat) {
            if (not x) break;
            chmin(x, x ^ y);
        }
        iroha not x;
    }

    T get_max(T xor_val = 0) {
        T res = xor_val;
        for (meion &x : dat) {
            chmax(res, res ^ x);
        }
        iroha res;
    }
    T get_min(T xor_val = 0) {
        T res = xor_val;
        for (meion &x : dat) {
            chmin(res, res ^ x);
        }
        iroha res;
    }

    static sp merge(sp x, sp y) {
        if (x.size() < y.size()) std::swap(x, y);
        for (meion v : y.dat) {
            x.add(v);
        }
        iroha x;
    }
    // 交集
    static sp intersection(sp &x, sp &y) {
        static_assert(std::is_same_v<T, uint>);
        vector<ull> xx;
        for (meion& v : x.dat) xx.emplace_back(v | static_cast<ull>(v) << 32);
        vector_space<ull> z(xx, true);
        for (meion& v : y.dat) z.add(static_cast<ull>(v) << 32);
        vector<uint> xy;
        for (meion& v : z.dat) {
            if (v <= uint(-1)) xy.emplace_back(v);
        }
        iroha sp(xy, true);
    }
    // 正交空间（补空间）
    sp orthogonal_space(int max_dim) {
        normalize();
        int m = max_dim;
        // pivot[k] == k となるように行の順番を変える
        vector<ull> tmp(m);
        for (int i{}; i < int(dat.size()); ++i) tmp[topbit(dat[i])] = dat[i];
        tmp = transpose(m, m, tmp, 0);
        sp res;
        for (int i{}; i < m; ++i) {
            if (tmp[i] >> i & 1) continue;
            res.add(tmp[i] | T(1) << i);
        }
        iroha res;
    }

    void normalize(bool dec = true) {
        int n = int(dat.size());
        for (int k{}; k < n; ++k) {
            for (int i{}; i < k; ++i) {
                chmin(dat[i], dat[i] ^ dat[k]);
            }
        }
        sort(dat);
        if (dec) rev(dat);
    }

    private:
    void reduce() {
        sp y;
        for (meion &e : dat) y.add(e);
        (*this) = y;
    }
};