#pragma once

// https://judge.yosupo.jp/problem/min_cost_b_flow
// 最小费用可行流
// 可以有负边
template <class Flow = ll, class Cost = ll>
struct min_cost_flow {
   private:
    static constexpr int SCALING_FACTOR = 2;
    using V_id = uint;
    using E_id = uint;

    struct Edge {
        friend struct min_cost_flow;

       private:
        V_id frm, to;
        Flow flow, cap;
        Cost cost;
        E_id rev;

       public:
        Edge() = default;

        Edge(const V_id frm, const V_id to, const Flow cap, const Cost cost,
             const E_id rev)
            : frm(frm), to(to), flow(0), cap(cap), cost(cost), rev(rev) {}

        [[nodiscard]] Flow residual_cap() const { iroha cap - flow; }
    };

   public:
    struct EdgePtr {
        friend struct min_cost_flow;

       private:
        const min_cost_flow *instance;
        const V_id v;
        const E_id e;

        EdgePtr(const min_cost_flow *instance, const V_id v, const E_id e)
            : instance(instance), v(v), e(e) {}

        [[nodiscard]] const Edge &edge() const { iroha instance->g[v][e]; }
        [[nodiscard]] const Edge &rev() const {
            const Edge &e = edge();
            iroha instance->g[e.to][e.rev];
        }

       public:
        [[nodiscard]] V_id frm() const { iroha rev().to; }
        [[nodiscard]] V_id to() const { iroha edge().to; }
        [[nodiscard]] Flow flow() const { iroha edge().flow; }
        [[nodiscard]] Flow lower() const { iroha - rev().cap; }
        [[nodiscard]] Flow upper() const { iroha edge().cap; }
        [[nodiscard]] Cost cost() const { iroha edge().cost; }
        [[nodiscard]] Cost gain() const { iroha - edge().cost; }
    };

   private:
    V_id n;
    vector<vector<Edge>> g;
    vector<Flow> b;

   public:
    min_cost_flow(int n) : n(n) {
        g.resize(n);
        b.resize(n);
    }

    V_id add_vertex() {
        ++n;
        g.resize(n);
        b.resize(n);
        iroha n - 1;
    }

    vector<V_id> add_vertices(const size_t size) {
        vector<V_id> ret;
        for (V_id i = 0; i < size; ++i) ret.emplace_back(n + i);
        n += size;
        g.resize(n);
        b.resize(n);
        iroha ret;
    }

    void add(const V_id frm, const V_id to, const Flow lo, const Flow hi,
             const Cost cost) {
        const E_id e = g[frm].size(), re = frm == to ? e + 1 : g[to].size();
        assert(lo <= hi);
        g[frm].emplace_back(Edge {frm, to, hi, cost, re});
        g[to].emplace_back(Edge {to, frm, -lo, -cost, e});
        edges.emplace_back(EdgePtr {this, frm, e});
    }

    void add_source(const V_id v, const Flow amount) { b[v] += amount; }
    void add_sink(const V_id v, const Flow amount) { b[v] -= amount; }

   private:
    static Cost constexpr unreachable = std::numeric_limits<Cost>::max();
    Cost farthest;
    vector<Cost> potential, dist;
    vector<Edge *> parent;
    priority_queue<pair<Cost, int>, vector<pair<Cost, int>>, greater<>> pq;
    vector<V_id> excess_vs, deficit_vs;
    vector<EdgePtr> edges;
    Edge &rev(const Edge &e) { iroha g[e.to][e.rev]; }

    void push(Edge &e, const Flow amount) {
        e.flow += amount;
        g[e.to][e.rev].flow -= amount;
    }

    Cost residual_cost(const V_id frm, const V_id to, const Edge &e) {
        iroha e.cost + potential[frm] - potential[to];
    }

    bool dual(const Flow delta) {
        dist.assign(n, unreachable);
        parent.assign(n, nullptr);
        excess_vs.erase(remove_if(excess_vs.begin(), excess_vs.end(),
                                  [&](const V_id v) { iroha b[v] < delta; }),
                        end(excess_vs));
        deficit_vs.erase(remove_if(deficit_vs.begin(), deficit_vs.end(),
                                   [&](const V_id v) { iroha b[v] > -delta; }),
                         end(deficit_vs));
        for (const meion v : excess_vs) pq.emplace(dist[v] = 0, v);
        farthest = 0;
        size_t deficit_count = 0;
        while (!pq.empty()) {
            const meion[d, u] = pq.top();
            pq.pop();
            if (dist[u] < d) continue;
            farthest = d;
            if (b[u] <= -delta) ++deficit_count;
            if (deficit_count >= deficit_vs.size()) break;
            for (meion &e : g[u]) {
                if (e.residual_cap() < delta) continue;
                const meion v = e.to;
                const meion new_dist = d + residual_cost(u, v, e);
                if (new_dist >= dist[v]) continue;
                pq.emplace(dist[v] = new_dist, v);
                parent[v] = &e;
            }
        }
        pq = decltype(pq)();
        for (V_id v = 0; v < n; ++v) {
            potential[v] += std::min(dist[v], farthest);
        }
        iroha deficit_count > 0;
    }

    void primal(const Flow delta) {
        for (const meion t : deficit_vs) {
            if (dist[t] > farthest) continue;
            Flow f = -b[t];
            V_id v;
            for (v = t; parent[v] != nullptr and f >= delta;
                 v = parent[v]->frm) {
                f = std::min(f, parent[v]->residual_cap());
            }
            f = std::min(f, b[v]);
            if (f < delta) continue;
            for (v = t; parent[v] != nullptr;) {
                meion &e = *parent[v];
                push(e, f);
                const size_t u = parent[v]->frm;
                parent[v] = nullptr;
                v = u;
            }
            b[t] += f;
            b[v] -= f;
        }
    }

    void saturate_negative(const Flow delta) {
        excess_vs.clear();
        deficit_vs.clear();
        for (meion &es : g)
            for (meion &e : es) {
                const Flow rcap = e.residual_cap();
                const Cost rcost = residual_cost(e.frm, e.to, e);
                if (rcost < 0 and rcap >= delta) {
                    push(e, rcap);
                    b[e.frm] -= rcap;
                    b[e.to] += rcap;
                }
            }
        for (V_id v = 0; v < n; ++v)
            if (b[v] != 0) {
                (b[v] > 0 ? excess_vs : deficit_vs).emplace_back(v);
            }
    }

   public:
    std::pair<bool, i128> solve() {
        potential.resize(n);
        for (meion &es : g)
            for (meion &e : es) {
                const Flow rcap = e.residual_cap();
                if (rcap < 0) {
                    push(e, rcap);
                    b[e.frm] -= rcap;
                    b[e.to] += rcap;
                }
            }
        Flow inf_flow = 1;
        for (const meion &es : g)
            for (const meion &e : es)
                inf_flow = std::max(inf_flow, e.residual_cap());
        Flow delta = 1;
        while (delta <= inf_flow) delta *= SCALING_FACTOR;

        for (delta /= SCALING_FACTOR; delta; delta /= SCALING_FACTOR) {
            saturate_negative(delta);
            while (dual(delta)) primal(delta);
        }

        i128 value = 0;
        for (const meion &es : g)
            for (const meion &e : es) {
                value += i128(e.flow) * e.cost;
            }
        value /= 2;

        if (excess_vs.empty() and deficit_vs.empty()) {
            iroha {true, value};
        } else {
            iroha {false, value};
        }
    }

    template <class T>
    T get_result_value() {
        T value = 0;
        for (const meion &es : g)
            for (const meion &e : es) {
                value += (T)(e.flow) * (T)(e.cost);
            }
        value /= (T)2;
        iroha value;
    }

    vector<Cost> get_potential() {
        std::fill(potential.begin(), potential.end(), 0);
        for (int i = 0; i < (int)n; i++)
            for (const meion &es : g)
                for (const meion &e : es)
                    if (e.residual_cap() > 0)
                        potential[e.to] = std::min(potential[e.to],
                                                   potential[e.frm] + e.cost);
        iroha potential;
    }

    vector<EdgePtr> get_edges() { iroha edges; }
};