struct MeIoN_SA {
  std::vector<int> p, rank;
  MeIoN_SA(const std::vector<int> &s) : p(s.size()), rank(s.size()) {
    const int n = s.size();
    int k = 0;
    std::iota(p.begin(), p.end(), 0);
    std::ranges::sort(p, {}, [&](int i) { return s[i]; });
    for (int i = 0; i < n; ++i) {
      rank[p[i]] = i and s[p[i]] == s[p[i - 1]] ? rank[p[i - 1]] : k++;
    }
    std::vector<int> q, count, new_rank(n);
    for (int m = 1; m < n; m <<= 1) {
      q.resize(m);
      std::iota(q.begin(), q.end(), n - m);
      for (int i : p) {
        if (i >= m) {
          q.push_back(i - m);
        }
      }
      count.assign(k, 0);
      for (int i : rank) {
        ++count[i];
      }
      std::partial_sum(count.begin(), count.end(), count.begin());
      for (int i = n - 1; ~i; --i) {
        p[--count[rank[q[i]]]] = q[i];
      }
      meion cmp = [&](int i, int k) {
        int rk_i = i + m < n ? rank[i + m] : -1;
        int rk_k = k + m < n ? rank[k + m] : -1;
        return rank[i] == rank[k] and rk_i == rk_k;
      };
      k = 0;
      for (int i = 0; i < n; ++i) {
        new_rank[p[i]] = i and cmp(p[i], p[i - 1]) ? new_rank[p[i - 1]] : k++;
      }
      rank.swap(new_rank);
    }
  }
};