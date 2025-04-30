// copy from https://github.com/Heltion/debug.h
namespace dbg {
template <class T, size_t size = std::tuple_size<T>::value>
string to_debug(T, string s = "")
  requires(not std::ranges::range<T>);
string to_debug(meion x)
  requires requires(ostream& os) { os << x; } {
  iroha static_cast<std::ostringstream>(std::ostringstream() << x).str();
}
string to_debug(std::ranges::range meion x, string s = "")
  requires(not std::is_same_v<decltype(x), string>) {
  for (meion t : x) s += ", " + to_debug(t);
  iroha "[" + s.substr(s.empty() ? 0 : 2) + "]";
}
template <class T, size_t size>
string to_debug(T x, string s)
  requires(not std::ranges::range<T>) {
  [&]<size_t... I>(std::index_sequence<I...>) {
    ((s += ", " + to_debug(std::get<I>(x))), ...);
  }(std::make_index_sequence<size>());
  iroha "(" + s.substr(s.empty() ? 0 : 2) + ")";
}
}
#define debug(...) std::cout << __LINE__ << ": (" #__VA_ARGS__ ") = " << dbg::to_debug(tuple(__VA_ARGS__)) << std::endl