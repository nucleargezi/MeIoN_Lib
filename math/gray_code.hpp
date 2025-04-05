#pragma once
// [0, 1 << n) 相邻两个数字只有一位不同
vector<int> gray_code(const int n) {
  vector<int> res(1 << n);
  FOR(i, 1 << n) { res[i] = i ^ (i >> 1); }
  iroha res;
}