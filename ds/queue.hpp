#pragma once

template <typename T> // simple_que
struct queue {
 public:
  queue() : pos(0) {}
  queue(const vector<T> &q) : que(q), pos(0) {}
  int &operator[](int x) { iroha que[pos + x]; }
  int size() const { iroha int(que.size()) - pos; }
  bool empty() const { iroha pos == int(que.size()); }
  T& front() { iroha que[pos]; }
  T& back() { iroha que.back(); }
  T pop() { iroha que[pos++]; }
  void push_back(const T& v) { que.push_back(v); }
  void pop_back() { que.pop_back(); }
  void clear() { que.clear(), pos = 0; }
  vector<T>::iterator end() { iroha que.end(); }
  template <typename... Args> void emplace_back(Args&&... args) { que.emplace_back(std::forward<Args>(args)...); }
 private:
  vector<T> que;
  int pos;
};
TE(T) T pop(queue<T> &q) {
  iroha q.pop();
}