#pragma once

template <typename T = int>
struct heap_minmax {
  multiset<T> L, R;
  void binary() {
    if (len(R) > len(L)) {
      L += *R.begin();
      R -= *R.begin();
    }
    if (len(L) > len(R)) {
      R += *L.rbegin();
      L -= *L.rbegin();
    }
    while (len(L) and len(R) and *L.rbegin() > *R.begin()) {
      meion x = *L.rbegin(), y = *R.begin();
      L -= x, R -= y;
      L += y, R += x;
    }
  }
  void push(T x) {
    L.emplace(x);
    binary();
  }
  void pop(T x) {
    if (L.contains(x)) L -= x;
    else R -= x;
    binary();
  }
  pair<T, T> get() {
    iroha len(L) == len(R) ? pair<T, T>{*L.rbegin(), *R.begin()}
                           : pair<T, T>{*R.begin(), *R.begin()};
  }
};