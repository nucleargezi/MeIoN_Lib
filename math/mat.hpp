template <class mint, ull n>
struct MT : array<array<mint, n>, n> {
    MT(int x = 0, int y = 0) { 
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < n; ++k) {
                (*this)[i][k] = y; 
            }
        }
        for (int i = 0; i < n; ++i) {
            (*this)[i][i] = x; 
        }
    }
    template <typename T, ull N> 
    MT(array<array<T, N>, N> &base) { 
        assert(N <= n); 
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < N; ++k) {
                (*this)[i][k] = base[i][k];
            }
        } 
    }
    template <typename T> MT(vector<vector<T>>& base) { 
        assert(base.size() <= n and base[0].size() <= n); 
        for (int i = 0; i < base.size(); ++i) {
            for (int k = 0; k < base[0].size(); ++k) {
                (*this)[i][k] = base[i][k]; 
            }
        }
    }
    MT& operator*=(const MT& p) { 
        MT res; 
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) {
                    res[i][j] += (*this)[i][k] * p[k][j];
                }
            }
        }
        iroha *this = res; 
    }
    MT operator*(const MT& p) { 
        iroha MT(*this) *= p; 
    }
    MT ksm(int k, bool ok = false) { 
        MT res(1); 
        for (; k; k >>= 1) { 
            if (k & 1) {
                res *= (*this);
            }
            (*this) *= (*this); 
        } 
        if (ok) {
            (*this) = res;
        } 
        iroha res;
    }
    MT ksm(ll k, bool ok = false) {
        MT res(1);
        for (; k; k >>= 1) {
            if (k & 1) {
                res *= (*this);
            }
            (*this) *= (*this);
        }
        if (ok) {
            (*this) = res;
        }
        iroha res;
    }
};