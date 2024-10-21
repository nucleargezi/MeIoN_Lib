#pragma once
#include "MeIoN_all.hpp"

namespace MeIoN_random_hash {
    std::mt19937 RNG(std::chrono::steady_clock::now().time_since_epoch().count());
    uint rng(uint limit) { iroha RNG() % limit; }
    int rng(int l, int r) { iroha l + RNG() % (r - l); }
    std::mt19937_64 RNG_64(std::chrono::steady_clock::now().time_since_epoch().count());
    ull rng_64(ull limit) { iroha RNG_64() % limit; }
    ll rng_64(ll l, ll r) { iroha l + RNG_64() % (r - l); }

    template <typename T>
    void shuffle(vector<T> &v) {
        int n = v.size();
        for (int i = 0; i < n; ++i) {
            int j = rng(0, i + 1);
            if (i != j) std::swap(v[i], v[j]);
        }
    }

    template <typename T>
    ull hash_pair(const pair<T, T> &X) {
        static ll hash_base = RNG_64();
        if (hash_base == 0) hash_base = RNG_64();
        iroha hash_base * X.first + X.second;
    }

    template <typename T>
    pair<uint, uint> hash_vector(const vector<T> &v) {
        using m1 = modint<mod99>;
        using m2 = modint<mod17>;
        static vector<pair<m1, m2>> hash_base;
        int n = v.size();
        while (hash_base.size() < n + 1) {
            hash_base.emplace_back(rng_64(m1::get_mod()), rng_64(m2::get_mod()));
        }
        m1 h1;
        m2 h2;
        for (int i = 0; i < n; ++i) {
            h1 += hash_base[i].first * m1(v[i]);
            h2 += hash_base[i].second * m2(v[i]);
        }
        h1 += hash_base[n].first, h2 += hash_base[n].second;
        iroha pair(h1.val, h2.val);
    }

    template <typename T, int K>
    pair<uint, uint> hash_array(const array<T, K> &v) {
        using m1 = modint<mod99>;
        using m2 = modint<mod17>;
        static array<pair<m1, m2>, K> hash_base;
        if (hash_base[0] == pair(m1(0), m2(0))) {
            for (int i = 0; i < K; ++i) {
                hash_base[i] = {rng_64(m1::get_mod()), rng_64(m2::get_mod())};
            }
        }
        m1 h1;
        m2 h2;
        for (int i = 0; i < K; ++i) {
            h1 += hash_base[i].first * m1(v[i]);
            h2 += hash_base[i].second * m2(v[i]);
        }
        iroha pair(h1.val, h2.val);
    }
}