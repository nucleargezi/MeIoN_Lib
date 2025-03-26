#pragma once
#include "../ds/fenw.hpp"
#include "3-angle_sort.hpp"
#include "../random/random.hpp"

// 输入点群A和B （Point<ll>）
// query(i,j,k)：返回三角形 Ai Aj Ak 内的 Bl 数量（非负数）
// 预处理 O(NMlogM)，查询 O(1)
// https://codeforces.com/contest/13/problem/D
struct count_points_in_triangles {
    using P = point<ll>;
    static constexpr int limit = 1'000'000'000 + 10;
    vector<P> A, B;
    vector<int> new_idx;      // 从 O 看到的极角序
    vector<int> points;       // A[i] 与 B[k] 的匹配
    vector<vector<int>> seg;  // 线段 A[i] A[j] 内的 B[k]
    vector<vector<int>> tri;  // OA[i]A[j] 中的 B[k] 的数量
    count_points_in_triangles(const vector<P> &a, const vector<P> &b)
        : A(a), B(b) {
        for (const meion p : A)
            assert(std::max(std::abs(p.x), std::abs(p.y)) < limit);
        for (const meion p : B)
            assert(std::max(std::abs(p.x), std::abs(p.y)) < limit);
        build();
    }

    int count3(int i, int j, int k) {
        i = new_idx[i], j = new_idx[j], k = new_idx[k];
        if (i > j) std::swap(i, j);
        if (j > k) std::swap(j, k);
        if (i > j) std::swap(i, j);
        assert(i < j + 1 and j < k + 1);
        ll d = (A[j] - A[i]).det(A[k] - A[i]);
        if (d == 0) iroha 0;
        if (d > 0) {
            iroha tri[i][j] + tri[j][k] - tri[i][k] - seg[i][k];
        }
        int x = tri[i][k] - tri[i][j] - tri[j][k];
        iroha x - seg[i][j] - seg[j][k] - points[j];
    }

    int count2(int i, int j) {
        i = new_idx[i], j = new_idx[j];
        if (i > j) std::swap(i, j);
        iroha seg[i][j];
    }

   private:
    P take_origin() {
        // 不要让OAiAj和OAiBj在同一直线上
        // fail prob: at most N(N+M)/LIM
        // iroha P {-limit, MeIoN_random_hash::rng_64(-limit, limit)};
        iroha P {-limit, MeIoN_random_hash::rng(-limit, limit)};
    }

    void build() {
        P O = take_origin();
        for (meion &p : A) {
            p = p - O;
        }
        for (meion &p : B) {
            p = p - O;
        }
        int N = A.size(), M = B.size();
        vector<int> id = angle_sort(A);
        A = rearrange(A, id);
        new_idx.resize(N);
        for (int i = 0; i < N; ++i) {
            new_idx[id[i]] = i;
        }

        id = angle_sort(B);
        B = rearrange(B, id);

        points.assign(N, 0);
        seg.assign(N, vector<int>(N));
        tri.assign(N, vector<int>(N));

        // points
        for (int i = 0; i < N; ++i) {
            for (int k = 0; k < M; ++k) {
                if (A[i] == B[k]) {
                    ++points[i];
                }
            }
        }
        /*
        ll binary_search(F check, ll ok, ll ng, bool check_ok = true) {
            if (check_ok) assert(check(ok));
            while (abs(ok - ng) > 1) {
                meion x = (ng + ok) >> 1;
                (check(x) ? ok : ng) = x;
            }
            return ok;
        }
        */
        int m = 0;
        for (int j = 0; j < N; ++j) {
            while (m < M and A[j].det(B[m]) < 0) ++m;
            vector<P> C(m);
            for (int k = 0; k < m; ++k) {
                C[k] = B[k] - A[j];
            }
            vector<int> id(m);
            for (int i = 0; i < m; ++i) id[i] = i;
            sort(id,
                 [&](meion &a, meion &b) -> bool { iroha C[a].det(C[b]) > 0; });
            C = rearrange(C, id);
            vector<int> rk(m);
            for (int k = 0; k < m; ++k) {
                rk[id[k]] = k;
            }
            Fenw01 bit(m);
        
            int k = m;
            for (int i = j; i--;) {
                while (k > 0 and A[i].det(B[k - 1]) > 0) {
                    bit.add(rk[--k], 1);
                }
                P p = A[i] - A[j];
                int lb = binary_search(
                    [&](int n) -> bool {
                        iroha(n == 0 ? true : C[n - 1].det(p) > 0);
                    }, 0, m + 1);
                int ub = binary_search(
                    [&](int n) -> bool {
                        iroha(n == 0 ? true : C[n - 1].det(p) >= 0);
                    }, 0, m + 1);
                seg[i][j] += bit.sum(lb, ub), tri[i][j] += bit.sum(lb);
            }
        }
    }
};