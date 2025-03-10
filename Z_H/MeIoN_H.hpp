#pragma once
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <random>
#include <ranges>
#include <set>
#include <string>
#include <tuple>
#include <utility>

using   std::array, std::bitset, std::deque, std::greater, std::less, std::map, 
        std::multiset, std::pair, std::priority_queue, std::set, 
        std::string, std::vector, std::tuple, std::function;

using NAME = void;       using uint = unsigned;   using ll = long long;      using ull = unsigned long long;     
using ld = long double;  using i128 = __int128;   using u128 = __uint128_t;  using f128 = __float128;

#define meion     auto
#define iroha     return
#define ST_T      auto  start = std::chrono::high_resolution_clock::now()
#define OT_T      auto  end = std::chrono::high_resolution_clock::now(); std::chrono::duration<double> elapsed = end - start; std::cout << "Elapsed time: " << elapsed.count() << "s\n"
#define TRUE(...)   (__VA_ARGS__, true)