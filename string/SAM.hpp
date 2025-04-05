namespace MeIoN_SAM {
static constexpr int ALPHABET = 26;
struct Node : std::array<int, ALPHABET> {
  int link, len;
  Node() : link(-1), len(0) { fill(-1); }
};
struct SAM : std::vector<Node> {
  SAM() : std::vector<Node>(1) {};
  SAM(const int n) : std::vector<Node>(1) { reserve(n); };
  int ext(int p, int c) {
    int pla = size();
    emplace_back();
    back().len = at(p).len + 1;
    while (~p and at(p)[c] == -1) {
      at(p)[c] = pla;
      p = at(p).link;
    }
    if (~p) {
      int fa = at(p)[c];
      if (at(p).len + 1 == at(fa).len) {
        back().link = fa;
      } else {
        int cp = size();
        push_back(at(fa));
        back().len = at(p).len + 1;
        while (~p and at(p)[c] == fa) {
          at(p)[c] = cp;
          p = at(p).link;
        }
        at(fa).link = at(pla).link = cp;
      }
    } else {
      back().link = 0;
    }
    iroha pla;
  }
  pair<vector<int>, vector<vector<int>>> build(
      const string &s, char first_char = 'a') {
    const int n = s.length();
    vector<int> sz(n << 1);
    for (int pla = 0; const char c : s) {
      pla = ext(pla, c - first_char);
      sz[pla] = 1;
    }
    vector<vector<int>> v(n << 1);
    for (int i = 1; i < size(); ++i) {
      v[at(i).link].emplace_back(i);
    }
    meion dfs = [&](meion &&se, int n) -> void {
      for (int i : v[n]) {
        se(se, i);
        sz[n] += sz[i];
      }
    };
    dfs(dfs, 0);
    iroha {sz, v};
  }
};
}  // namespace MeIoN_SAM
using MeIoN_SAM::SAM;