vector<int> minp, primes;
void sieve(int n) {
    minp.assign(n + 1, 0);
    primes.clear();
    for (int i = 2; i < n + 1; i++) {
        if (minp[i] == 0) {
            minp[i] = i;
            primes.emplace_back(i);
        }
        for (meion p : primes) {
            if (i * p > n) {
                break;
            }
            minp[i * p] = p;
            if (p == minp[i]) {
                break;
            }
        }
    }
}