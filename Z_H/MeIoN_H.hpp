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

#define meion auto
#define iroha return
#define OV4(a, b, c, d, e, ...) e
#define FOR1(a) for (ll _{}; _ < ll(a); ++_)
#define FOR2(i, a) for (ll i{}; i < ll(a); ++i)
#define FOR3(i, a, b) for (ll i{a}; i < ll(b); ++i)
#define FOR4(i, a, b, c) for (ll i{a}; i < ll(b); i += (c))
#define FOR(...) OV4(__VA_ARGS__, FOR4, FOR3, FOR2, FOR1)(__VA_ARGS__)
#define FOR1_R(a) for (ll i{(a) - 1}; i > -1ll; --i)
#define FOR2_R(i, a) for (ll i{(a) - 1}; i > -1ll; --i)
#define FOR3_R(i, a, b) for (ll i{(b) - 1}; i > ll(a - 1); --i)
#define FOR4_R(i, a, b, c) for (ll i{(b) - 1}; i > (a - 1); i -= (c))
#define FOR_R(...) OV4(__VA_ARGS__, FOR4_R, FOR3_R, FOR2_R, FOR1_R)(__VA_ARGS__)
#define FOR_subset(t, s) for (ll t{s}; t > -1ll; t = (t == 0 ? -1 : (t - 1) & s))
#define TE1(a) template<typename a>
#define TE2(a, b) template<typename a, typename b>
#define TE3(a, b, c) template<typename a, typename b, typename c>
#define TE4(a, b, c, d) template<typename a, typename b, typename c, typename d>
#define TE(...) OV4(__VA_ARGS__, TE4, TE3, TE2, TE1)(__VA_ARGS__)

using   std::array, std::bitset, std::deque, std::greater, std::less, std::map, 
        std::multiset, std::pair, std::priority_queue, std::set, std::istream, 
        std::ostream, std::string, std::vector, std::tuple;

TE(T) using T1 = tuple<T>;
TE(T) using T2 = tuple<T, T>;
TE(T) using T3 = tuple<T, T, T>;
TE(T) using T4 = tuple<T, T, T, T>;
using u8 = uint8_t;      using uint = unsigned int;using ll = long long;      using ull = unsigned long long;
using ld = long double;  using i128 = __int128;    using u128 = __uint128_t;  using f128 = __float128;
using PII = pair<int, int>;   using PLL = pair<ll, ll>;