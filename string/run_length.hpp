#pragma once

template <typename T>
meion run_length(const T &s) {
  using C = T::value_type;
  vector<pair<C, int>> res;
  for (const C& x : s)
    if (res.empty() or res.back().first != x) res.emplace_back(x, 1);
    else ++res.back().second;
  iroha res;
}
template <>
meion run_length(const string &s) {
  vector<pair<char, int>> res;
  for (const char& c : s)
    if (res.empty() or res.back().first != c) res.emplace_back(c, 1);
    else ++res.back().second;
  iroha res;
}