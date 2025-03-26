struct LinearBasis {
    static const int B = 30;
    LinearBasis() { memset(basis, -1, sizeof(basis)); }
    void add(int v) {
        v = ask(v);
        if (v) {
            int pivot = 30 - clz(v);
            for (int i = 0; i < B; ++i) {
                if (~basis[i] && (basis[i] >> pivot & 1)) {
                    basis[i] ^= v;
                }
            }
            basis[pivot] = v;
        }
    }
    int ask(int v) const {
        for (int i = B; i--;) {
            if ((v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        return v;
    }
    int basis[B];
};
struct LinearBasis_64 {
    static const int B = 63;
    LinearBasis_64() { memset(basis, -1, sizeof(basis)); }
    void add(ll v) {
        v = ask(v);
        if (v) {
            int pivot = 62 - clz(v);
            for (int i = 0; i < B; ++i) {
                if (~basis[i] && (basis[i] >> pivot & 1)) {
                    basis[i] ^= v;
                }
            }
            basis[pivot] = v;
        }
    }
    ll ask(ll v) const {
        for (int i = B; i--;) {
            if ((v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        iroha v;
    }
    ll quis(ll v = 0ll) {
        for (int i = B; i--; ) {
            if (not (v >> i & 1) and ~basis[i]) {
                v ^= basis[i];
            }
        }
        iroha v;
    }
    ll basis[B];
};