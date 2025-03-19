// 似乎在某些情况下会 死循环 / 随机越界...
template <typename Cap = int>
struct HLPP {
   public:
    struct edge {
        int to, rev;
        Cap flow;
        edge(int to, Cap flow, int rev) : to(to), flow(flow), rev(rev) {}
        edge() {}
    };
    int n, m, s, t, maxh, maxgaph, workcnt;
    std::vector<std::vector<edge>> vec;
    std::vector<Cap> ov;
    std::vector<int> h, cur, ovlist, ovnex, gap, gapprv, gapnex;
    HLPP(int n, int m, int s, int t)
        : n(n),
          m(m),
          s(s),
          t(t),
          maxh(0),
          maxgaph(0),
          workcnt(0),
          vec(n),
          ov(n),
          h(n),
          cur(n),
          ovlist((n), -1),
          ovnex((n), -1),
          gap((n), -1),
          gapprv((n), -1),
          gapnex((n), -1) {}
    HLPP() {}
    void add(int u, int v, Cap f) {
        vec[u].push_back(edge(v, f, vec[v].size()));
        vec[v].push_back(edge(u, 0, vec[u].size() - 1));
    }

    Cap flow() {
        relabel();
        for (auto &e : vec[s]) {
            if (e.flow) {
                push_flow(s, e, e.flow);
                maxh = std::max(maxh, h[e.to]);
            }
        }
        for (; maxh; --maxh) {
            while (~ovlist[maxh]) {
                int x = ovlist[maxh];
                ovlist[maxh] = ovnex[x];
                update_dis(x);
                if (workcnt >= (n << 2)) {
                    relabel();
                }
            }
        }
        return ov[t];
    }

    void update_dis(int x) {
        int nh = n, sz = vec[x].size();
        for (int i = cur[x]; i < sz; ++i) {
            auto &e = vec[x][i];
            if (e.flow) {
                if (h[x] == h[e.to] + 1) {
                    push_flow(x, e, std::min(ov[x], e.flow));
                    if (!ov[x]) {
                        cur[x] = i;
                        return;
                    }
                } else {
                    nh = std::min(nh, h[e.to] + 1);
                }
            }
        }
        for (int i = 0; i < cur[x]; ++i) {
            auto &e = vec[x][i];
            if (e.flow) {
                nh = std::min(nh, h[e.to] + 1);
            }
        }
        cur[x] = 0;
        ++workcnt;
        if (~gapnex[gap[h[x]]]) {
            set_height(x, nh);
        } else {
            int oldh = h[x];
            for (int i = oldh; i <= maxgaph; ++i) {
                for (int j = gap[i]; ~j; j = gapnex[j]) {
                    h[j] = n;
                }
                gap[i] = -1;
            }
            maxgaph = oldh - 1;
        }
    }
    void set_height(int x, int newh) {
        if (~gapprv[x]) {
            if (gapprv[x] == x) {
                gapprv[gapnex[x]] = gapnex[x];
                gap[h[x]] = gapnex[x];
            } else {
                gapnex[gapprv[x]] = gapnex[x];
                if (~gapnex[x]) {
                    gapprv[gapnex[x]] = gapprv[x];
                }
            }
        }
        if ((h[x] = newh) >= n) {
            return;
        }
        maxgaph = std::max(maxgaph, newh);
        if (ov[x]) {
            maxh = std::max(maxh, h[x]);
            ovnex[x] = ovlist[h[x]];
            ovlist[h[x]] = x;
        }
        if (~(gapnex[x] = gap[h[x]])) {
            gapprv[gapnex[x]] = x;
        }
        gap[h[x]] = gapprv[x] = x;
    }
    void push_flow(int u, edge &e, Cap df) {
        if (!ov[e.to] && (e.to) != t) {
            ovnex[e.to] = ovlist[h[e.to]];
            ovlist[h[e.to]] = (e.to);
        }
        ov[u] -= df, ov[e.to] += df, e.flow -= df, vec[e.to][e.rev].flow += df;
    }
    void relabel() {
        workcnt = maxh = maxgaph = 0;
        std::fill(h.begin(), h.end(), n), h[t] = 0;
        std::fill(gap.begin(), gap.end(), -1);
        std::fill(gapprv.begin(), gapprv.end(), -1);
        std::fill(gapnex.begin(), gapnex.end(), -1);
        std::fill(ovlist.begin(), ovlist.end(), -1);
        std::fill(ovnex.begin(), ovnex.end(), -1);
        std::fill(cur.begin(), cur.end(), 0);
        std::queue<int> q;
        q.push(t);
        int x;
        while (!q.empty()) {
            x = q.front(), q.pop();
            for (auto &e : vec[x]) {
                if (h[e.to] == n && e.to != s && vec[e.to][e.rev].flow) {
                    set_height(e.to, h[x] + 1);
                    q.push(e.to);
                }
            }
        }
    }
};