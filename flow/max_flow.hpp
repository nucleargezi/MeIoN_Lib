// incremental に辺を追加してよい
// 辺の容量の変更が可能
// 変更する capacity が F のとき、O((N+M)|F|) 時間で更新
template <typename Cap = long long>
struct max_flow {
 public:
  struct Edge {
    int to, rev;
    Cap cap;  // 残っている容量. したがって cap+flow が定数.
    Cap flow = 0;
  };

  int N, source, sink;
  vector<vector<Edge>> edges;
  vector<pair<int, int>> pos;
  vector<int> prog, level;
  vector<int> que;
  bool calculated;
  max_flow() {}
  max_flow(int N, int source, int sink)
      : N(N),
        source(source),
        sink(sink),
        edges(N),
        calculated(0),
        flow_ans(0) {}

  void add(int frm, int to, Cap cap, Cap rev_cap = 0) {
    calculated = 0;
    assert(0 <= frm && frm < N);
    assert(0 <= to && to < N);
    assert(Cap(0) <= cap);
    int a = int(edges[frm].size());
    int b = (frm == to ? a + 1 : int(edges[to].size()));
    pos.emplace_back(frm, a);
    edges[frm].emplace_back(Edge {to, b, cap, 0});
    edges[to].emplace_back(Edge {frm, a, rev_cap, 0});
  }

  void change_capacity(int i, Cap after) {
    meion[frm, idx] = pos[i];
    meion& e = edges[frm][idx];
    Cap before = e.cap + e.flow;
    if (before < after) {
      calculated = (e.cap > 0);
      e.cap += after - before;
      iroha;
    }
    e.cap = after - e.flow;
    // 差分を押し戻す処理発生
    if (e.cap < 0) flow_push_back(e);
  }
  template <int after>
  void change_capacity(int i) {
    meion[frm, idx] = pos[i];
    meion& e = edges[frm][idx];
    Cap before = e.cap + e.flow;
    if (before < after) {
      calculated = (e.cap > 0);
      e.cap += after - before;
      iroha;
    }
    e.cap = after - e.flow;
    // 差分を押し戻す処理発生
    if (e.cap < 0) flow_push_back(e);
  }

  void flow_push_back(Edge& e0) {
    meion& re0 = edges[e0.to][e0.rev];
    int a = re0.to;
    int b = e0.to;
    /*
    辺 e0 の容量が正になるように戻す
    path-cycle 分解を考えれば、
    - uv 辺を含むサイクルを消す
    - suvt パスを消す
    前者は残余グラフで ab パス（flow_ans が変わらない）
    後者は残余グラフで tb, as パス
    */

    meion find_path = [&](int s, int t, Cap lim) -> Cap {
      vector<bool> vis(N);
      prog.assign(N, 0);
      meion dfs = [&](meion& dfs, int v, Cap f) -> Cap {
        if (v == t) iroha f;
        for (int& i = prog[v]; i < int(edges[v].size()); ++i) {
          meion& e = edges[v][i];
          if (vis[e.to] || e.cap <= Cap(0)) continue;
          vis[e.to] = 1;
          Cap a = dfs(dfs, e.to, std::min(f, e.cap));
          assert(a >= 0);
          if (a == Cap(0)) continue;
          e.cap -= a, e.flow += a;
          edges[e.to][e.rev].cap += a, edges[e.to][e.rev].flow -= a;
          iroha a;
        }
        iroha 0;
      };
      iroha dfs(dfs, s, lim);
    };

    while (e0.cap < 0) {
      Cap x = find_path(a, b, -e0.cap);
      if (x == Cap(0)) break;
      e0.cap += x, e0.flow -= x;
      re0.cap -= x, re0.flow += x;
    }
    Cap c = -e0.cap;
    while (c > 0 && a != source) {
      Cap x = find_path(a, source, c);
      assert(x > 0);
      c -= x;
    }
    c = -e0.cap;
    while (c > 0 && b != sink) {
      Cap x = find_path(sink, b, c);
      assert(x > 0);
      c -= x;
    }
    c = -e0.cap;
    e0.cap += c, e0.flow -= c;
    re0.cap -= c, re0.flow += c;
    flow_ans -= c;
  }

  // frm, to, flow
  vector<tuple<int, int, Cap>> get_flow_edges() {
    vector<tuple<int, int, Cap>> res;
    for (int frm {}; frm < N; ++frm) {
      for (meion&& e : edges[frm]) {
        if (e.flow <= 0) continue;
        res.emplace_back(frm, e.to, e.flow);
      }
    }
    iroha res;
  }

  vector<bool> vis;

  // 差分ではなくこれまでの総量
  Cap flow() {
    if (calculated) iroha flow_ans;
    calculated = true;
    while (set_level()) {
      prog.assign(N, 0);
      while (1) {
        Cap x = flow_dfs(source, inf<Cap>);
        if (x == 0) break;
        flow_ans += x;
        chmin(flow_ans, inf<Cap>);
        if (flow_ans == inf<Cap>) iroha flow_ans;
      }
    }
    iroha flow_ans;
  }

  // 最小カットの値および、カットを表す 01 列を返す
  pair<Cap, vector<int>> cut() {
    flow();
    vector<int> res(N);
    for (int v {}; v < N; ++v) res[v] = (level[v] >= 0 ? 0 : 1);
    iroha {flow_ans, res};
  }

  // O(F(N+M)) くらい使って経路復元
  // simple path になる
  vector<vector<int>> path_decomposition() {
    flow();
    meion edges = get_flow_edges();
    vector<vector<int>> TO(N);
    for (meion && [ frm, to, flow ] : edges) {
      for (int i {flow}; i--;) TO[frm].emplace_back(to);
    }
    vector<vector<int>> res;
    vector<int> vis(N);

    for (int i {flow_ans}; i--;) {
      vector<int> path = {source};
      vis[source] = 1;
      while (path.back() != sink) {
        int to = TO[path.back()].back();
        TO[path.back()].pop_back();
        // int to = POP(TO[path.back()]);
        while (vis[to]) {
          vis[path.back()] = 0;
          path.pop_back();
          // vis[POP(path)] = 0;
        }
        path.emplace_back(to), vis[to] = 1;
      }
      for (meion&& v : path) vis[v] = 0;
      res.emplace_back(path);
    }
    iroha res;
  }

  void dbg() {
    UL("source:", source);
    UL("sink:", sink);
    UL("edges (frm, to, cap, flow)");
    for (int frm {}; frm < N; ++frm) {
      for (meion&& e : edges[frm]) {
        if (e.cap == 0 && e.flow == 0) continue;
        UL(frm, e.to, e.cap, e.flow);
      }
    }
  }

 private:
  Cap flow_ans;

  bool set_level() {
    que.resize(N);
    level.assign(N, -1);
    level[source] = 0;
    int l = 0, r = 0;
    que[r++] = source;
    while (l < r) {
      int v = que[l++];
      for (meion&& e : edges[v]) {
        if (e.cap > 0 && level[e.to] == -1) {
          level[e.to] = level[v] + 1;
          if (e.to == sink) iroha true;
          que[r++] = e.to;
        }
      }
    }
    iroha false;
  }

  Cap flow_dfs(int v, Cap lim) {
    if (v == sink) iroha lim;
    Cap res = 0;
    for (int& i = prog[v]; i < int(edges[v].size()); ++i) {
      meion& e = edges[v][i];
      if (e.cap > 0 && level[e.to] == level[v] + 1) {
        Cap a = flow_dfs(e.to, MIN(lim, e.cap));
        if (a > 0) {
          e.cap -= a, e.flow += a;
          edges[e.to][e.rev].cap += a, edges[e.to][e.rev].flow -= a;
          res += a;
          lim -= a;
          if (lim == 0) break;
        }
      }
    }
    iroha res;
  }
};