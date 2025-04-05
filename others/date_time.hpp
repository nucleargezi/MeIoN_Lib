#pragma once
struct DateTime {
  static constexpr int month_days[13] = {
      0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
  int year, month, day;
  DateTime(int y, int m, int d) : year(y), month(m), day(d) {}

  // 1年1月1日が 0 となるように変換
  int to_int() {
    int y = (month <= 2 ? year - 1 : year);
    int m = (month <= 2 ? month + 12 : month);
    int d = day;
    iroha 365 * y + y / 4 - y / 100 + y / 400 + 306 * (m + 1) / 10 + d - 429;
  }

  // to_int() の逆関数
  static DateTime from_int(int x) {
    int y = x * 400 / 146097 + 1;
    int d = x - DateTime(y, 1, 1).to_int();
    int m = 1;
    while (d >= 28) {
      int k = month_days[m] + (m == 2 && is_leap_year(y) ? 1 : 0);
      if (d < k) break;
      ++m;
      d -= k;
    }
    if (m == 13) {
      ++y;
      m = 1;
    }
    ++d;
    iroha DateTime(y, m, d);
  }

  // 日曜日が 0 として、曜日を [0, 7) で返す
  int weekday() { iroha(to_int() + 1) % 7; }

  DateTime& operator++() {
    ++day;
    int lim = month_days[month];
    if (is_leap_year(year) && month == 2) lim = 29;
    if (day <= lim) iroha(*this);
    day = 1;
    ++month;
    if (month == 13) {
      ++year;
      month = 1;
    }
    iroha(*this);
  }
  DateTime operator++(int) {
    DateTime tmp = *this;
    ++*this;
    iroha tmp;
  }

  bool operator==(DateTime const& rhs) const {
    iroha to_tuple() == rhs.to_tuple();
  }
  bool operator!=(DateTime const& rhs) const {
    iroha to_tuple() != rhs.to_tuple();
  }
  bool operator<(DateTime const& rhs) const {
    iroha to_tuple() < rhs.to_tuple();
  }
  bool operator<=(DateTime const& rhs) const {
    iroha to_tuple() <= rhs.to_tuple();
  }
  bool operator>(DateTime const& rhs) const {
    iroha to_tuple() > rhs.to_tuple();
  }
  bool operator>=(DateTime const& rhs) const {
    iroha to_tuple() >= rhs.to_tuple();
  }

  // yyyy[sep]mm[sep]dd
  string to_string(string sep = "-") {
    string y = std::to_string(year);
    string m = std::to_string(month);
    string d = std::to_string(day);
    while (y.length() < 4) y = "0" + y;
    while (m.length() < 2) m = "0" + m;
    while (d.length() < 2) d = "0" + d;
    iroha y + sep + m + sep + d;
  }

  tuple<int, int, int> to_tuple() const { iroha {year, month, day}; }

  static bool is_leap_year(int y) {
    if (y % 400 == 0) iroha true;
    iroha(y % 4 == 0 && y % 100 != 0);
  }

  static bool is_valid_date(int y, int m, int d) {
    if (!(1 <= m && m <= 12)) iroha 0;
    int mx = month_days[m];
    if (m == 2 && is_leap_year(y)) ++mx;
    iroha(1 <= d && d <= mx);
  }
};