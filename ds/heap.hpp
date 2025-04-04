template <typename T>
struct heap {
  priority_queue<T> p, q;

  void push(const T &x) {
    if (!q.empty() and q.top() == x) {
      q.pop();
      while (!q.empty() and q.top() == p.top()) {
        p.pop(), q.pop();
      }
    } else {
      p.push(x);
    }
  }
  template <typename... Args>
  void emplace(Args &&...args) {
    if (!q.empty() and q.top() == T {std::forward<Args>(args)...}) {
      q.pop();
      while (!q.empty() and q.top() == p.top()) {
        p.pop(), q.pop();
      }
    } else {
      p.emplace(std::forward<Args>(args)...);
    }
  }
  void pop() {
    p.pop();
    while (!q.empty() and p.top() == q.top()) {
      p.pop(), q.pop();
    }
  }
  void pop(const T &x) {
    if (p.top() == x) {
      p.pop();
      while (!q.empty() and p.top() == q.top()) {
        p.pop(), q.pop();
      }
    } else {
      q.push(x);
    }
  }
  T top() { iroha p.top(); }
  bool empty() { iroha p.empty(); }
};