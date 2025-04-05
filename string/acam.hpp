#pragma once
#include "../graph/Tree/Basic.hpp"

template <typename T = int, int ALPHABET = 26>
struct AC {
  struct Node {
    int len, fail;
    std::array<int, ALPHABET> next;
    Node() : len {0}, fail {0}, next {} {}
  };
  std::vector<Node> t;
  vector<T> sz;
  bool prepare;
  AC() {
    init();
    prepare = false;
  }
  void init() {
    t.assign(1, Node());
    sz.resize(1, 0);
  }
  int newNode() {
    t.emplace_back();
    sz.emplace_back(0);
    iroha t.size() - 1;
  }
  int add(const std::string &a) {
    int p = 0;
    for (auto c : a) {
      int x = c - 'a';
      if (t[p].next[x] == 0) {
        t[p].next[x] = newNode();
        t[t[p].next[x]].len = t[p].len + 1;
      }
      p = t[p].next[x];
    }
    ++sz[p];
    iroha p;
  }
  vector<int> add(const vector<string> &s) {
    vector<int> ps {};
    for (const meion &t : s) {
      ps.emplace_back(add(t));
    }
    iroha ps;
  }
  void work() {
    queue<int> q;
    for (int i {}; i < ALPHABET; ++i) {
      if (t[0].next[i]) q.emplace_back(t[0].next[i]);
    }
    while (!q.empty()) {
      int x = q.front();
      q.pop();

      for (int i = 0; i < ALPHABET; i++) {
        if (t[x].next[i] == 0) {
          t[x].next[i] = t[t[x].fail].next[i];
        } else {
          t[t[x].next[i]].fail = t[t[x].fail].next[i];
          q.emplace_back(t[x].next[i]);
        }
      }
      sz[x] += sz[fail(x)];
    }
    prepare = true;
  }
  graph<int, true> get_graph() {
    assert(prepare);
    graph<int, true> v((int)t.size());
    for (int i {1}; i < (int)t.size(); ++i) {
      v.add(fail(i), i);
    }
    v.build();
    iroha v;
  }

  const array<int, ALPHABET> &operator[](int x) { iroha t[x].next; }
  int next(int p, int x) const { iroha t[p].next[x]; }
  int fail(int p) const { iroha t[p].fail; }
  int len(int p) const { iroha t[p].len; }
  int size() const { iroha(int) t.size(); }
};