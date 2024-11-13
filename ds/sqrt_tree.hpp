#pragma once
template <typename Monoid>
struct sqrt_tree {   // nlog^2 预处理 O1查询区间信息 满足结合律
    using MX = Monoid;
    using X = typename MX::value_type;

    static constexpr int K = 3;
    static constexpr uint SZ[] = {8, 64, 4096};
    static constexpr uint MASK[] = {7, 63, 4095};

    int N;
    // 元となる静的な列
    vector<X> A;
    // 各階層に対して，ブロック先頭からある要素まで [s,i]
    // 各階層に対して，ある要素からブロック末尾まで [i,t]
    vector<vector<X>> PREF, SUFF;
    // 各階層に対して，あるブロックからあるブロックまで
    vector<vector<X>> BETWEEN;

    sqrt_tree() {}
    template <typename F>
    sqrt_tree(int n, F f) {
        build(n, f);
    }
    sqrt_tree(const vector<X>& v) {
        build(v.size(), [&](int i) -> X { return v[i]; });
    }

    template <typename F>
    void build(int n_, F f) {
        N = n_;
        assert(N <= (1 << 24));
        A.reserve(N);
        for (int i = 0; i < N; ++i) A.emplace_back(f(i));
        // まず prefix, suffix の構築
        PREF.assign(K, A), SUFF.assign(K, A);
        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < N; ++i) {
                if (i & MASK[k]) PREF[k][i] = MX::op(PREF[k][i - 1], A[i]);
            }
            for (int i = N; --i; ) {
                if (i & MASK[k]) SUFF[k][i - 1] = MX::op(A[i - 1], SUFF[k][i]);
            }
        }
        // between の構築
        BETWEEN.resize(K);
        for (int k = 0; k < K; ++k) {
            // n : 全体の小ブロックの個数
            auto get = [&](int i) -> X { return SUFF[k][SZ[k] * i]; };
            int n = N / SZ[k];
            int s = 0;
            for (int r = 0; r < n; ++r) {
                if (r % SZ[k] == 0) s = r;
                BETWEEN[k].emplace_back(get(r));
                for (int l = r - 1; l >= s; --l) {
                    BETWEEN[k].emplace_back(MX::op(get(l), BETWEEN[k].back()));
                }
            }
        }
    }

    static constexpr int BIT_TO_LAYER[] = {0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2, 2,
                                           3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3};

    X prod(int L, int R) {
        assert(0 <= L && L <= R && R <= N);
        if (L == R) return MX::unit();
        if (L + 1 == R) return A[L];
        --R;
        int k = BIT_TO_LAYER[topbit(L ^ R)];
        if (k == 0) {
            // 長さ SZ[0] のブロックにクエリが収まっている. 愚直に.
            X x = A[L];
            for (int i = L + 1; i < R + 1; ++i) x = MX::op(x, A[i]);
            return x;
        }
        --k;
        // 同じ長さ SZ[k+1] のブロック内にある. 違う SZ[k] ブロック内にある.
        uint a = L / SZ[k], b = R / SZ[k];
        assert(a < b);
        X &x1 = SUFF[k][L], &x2 = PREF[k][R];
        if (a + 1 == b) return MX::op(x1, x2);
        ++a, --b;
        // [a,b] 番目の SZ[k]-block の間を取得する
        // BETWEEN のどこにデータが置いてあるか調べる
        uint m = a / SZ[k];
        a &= MASK[k], b &= MASK[k];
        uint idx = m * (SZ[k] / 2) * (SZ[k] + 1);
        idx += (b + 1) * (b + 2) / 2 - 1 - a;
        return MX::op(x1, MX::op(BETWEEN[k][idx], x2));
    }
};