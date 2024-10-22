namespace RMQ {
    vector<int> lg(2);
    template <typename T> struct maxtable {
        vector<T> a;
        vector<vector<T>> st;
        static int x;
        maxtable(const vector<T> &b):a(b) {
            int n = a.size(), i, j, k, r;
            while (lg.size()<=n) lg.emplace_back(lg[lg.size() >> 1] + 1);
            st.assign(lg[n] + 1,vector<T>(n));
            std::iota(st[0].begin(), st[0].end(), 0);
            for (j = 1; j <= lg[n]; j++) {
                r = n - (1 << j);
                k = 1 << j - 1;
                for (i = 0; i <= r; i++) 
                    st[j][i] = a[st[j-1][i]] < a[st[j-1][i+k]] ? st[j-1][i+k] : st[j-1][i];
            }
        }
        T quis(int l, int r) const { // [l, r]
            assert(0 <= l and l <= r and r < a.size());
            int z = lg[r - l + 1];
            return std::max(a[st[z][l]], a[st[z][r - (1 << z) + 1]]);
        }
        int rmp(int l,int r) const {
            assert(0 <= l and l <= r and r < a.size());
            int z = lg[r-l+1];
            return a[st[z][l]] < a[st[z][r - (1 << z) + 1]] ? st[z][r - (1 << z) + 1] : st[z][l];
        }
    };
} using RMQ::maxtable;