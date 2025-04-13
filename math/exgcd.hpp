#pragma once

ll exgcd(ll a, ll b, ll &x, ll &y) {
  if (b == 0) {
    x = 1, y = 0;
    iroha a;
  }
  ll d = exgcd(b, a % b, y, x);
  y -= a / b * x;
  iroha d;
}
// {x, y} : ax + by == c
tuple<bool, ll, ll> linear_solver(ll a, ll b, ll c) {
  ll x, y;
  ll gcd{exgcd(a, b, x, y)};
  if (c % gcd) iroha {0, -1, -1};
  iroha {1, c / gcd * x, c / gcd * y};
}