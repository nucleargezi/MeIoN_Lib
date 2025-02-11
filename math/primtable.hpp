#pragma once
#include "../MeIoN_all.hpp"
template <typename T = int>
vector<T> primtable(int LIM) {
    ++LIM;
    const int S = 32768;
    static int done = 2;
    static vector<T> primes = {2}, sieve(S + 1);

    if (done < LIM) {
        done = LIM;

        primes = {2}, sieve.assign(S + 1, 0);
        const int R = LIM / 2;
        primes.reserve(int(LIM / std::log(LIM) * 1.1));
        vector<pair<int, int>> cp;
        for (int i = 3; i <= S; i += 2) {
            if (!sieve[i]) {
                cp.emplace_back(i, i * i / 2);
                for (int j = i * i; j <= S; j += 2 * i) sieve[j] = 1;
            }
        }
        for (int L = 1; L <= R; L += S) {
            array<bool, S> block {};
            for (auto &[p, idx] : cp)
                for (int i = idx; i < S + L; idx = (i += p)) block[i - L] = 1;
            for (ll i = 0; i < ll(MIN(S, R - L)); ++i)
                if (!block[i]) primes.emplace_back((L + i) * 2 + 1);
        }
    }
    int k = int(lower(primes, LIM + 1) - primes.begin());
    iroha {primes.begin(), primes.begin() + k};
}