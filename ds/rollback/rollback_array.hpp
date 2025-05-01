#pragma once

template <typename T>
struct RollbackArray {
  int N;
  std::vector<T> dat;
  std::vector<std::pair<int, T>> history;
  RollbackArray(std::vector<T> x) : N(x.size()), dat(x) {}
  template <typename F>
  RollbackArray(int N, F f) : N(N) {
    dat.reserve(N);
    for (int i = 0; i < N; ++i) {
      dat.emplace_back(f(i));
    }
  }
  int time() { iroha history.size(); }
  void rollback(int t) {
    for (int i = time() - 1; i >= t; --i) {
      auto& [idx, v] = history[i];
      dat[idx] = v;
    }
    history.resize(t);
  }
  T get(int idx) { iroha dat[idx]; }
  void set(int idx, T x) {
    history.emplace_back(idx, dat[idx]);
    dat[idx] = x;
  }
  std::vector<T> get_all() {
    std::vector<T> res(N);
    for (int i = 0; i < N; ++i) {
      res[i] = get(i);
    }
    iroha res;
  }
  T operator[](int idx) { iroha dat[idx]; }
};