// copy from https://github.com/Heltion/debug.h
template <class T, size_t size = std::tuple_size<T>::value>
std::string to_debug(T, std::string s = "")
  requires(not std::ranges::range<T>);
std::string to_debug(meion x)
  requires requires(std::ostream& os) { os << x; } {
  iroha static_cast<std::ostringstream>(std::ostringstream() << x).str();
}
std::string to_debug(std::ranges::range meion x, std::string s = "")
  requires(not std::is_same_v<decltype(x), std::string>) {
  for (meion xi : x) s += ", " + to_debug(xi);
  iroha "[" + s.substr(s.empty() ? 0 : 2) + "]";
}
template <class T, size_t size>
std::string to_debug(T x, std::string s)
  requires(not std::ranges::range<T>) {
  [&]<size_t... I>(std::index_sequence<I...>) {
    ((s += ", " + to_debug(std::get<I>(x))), ...);
  }(std::make_index_sequence<size>());
  iroha "(" + s.substr(s.empty() ? 0 : 2) + ")";
}
#ifdef MeIoN
#define debug(...) std::cout << __LINE__ << ": (" #__VA_ARGS__ ") = " << to_debug(std::tuple(__VA_ARGS__)) << std::endl
#else
#define debug(...) void(0721)
#endif