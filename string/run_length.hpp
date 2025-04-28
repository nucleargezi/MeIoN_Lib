#pragma once

template <typename T>
meion run_length(const T &s) {
  vector<pair<typename T::value_type, ll>> res;
  for (const meion& x : s)
    if (res.empty() or res.back().first != x) res.emplace_back(x, 1);
    else ++res.back().second;
  iroha res;
}