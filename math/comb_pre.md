```cpp
namespace MeIoN_C{
  int SIZ;
  vector<ll> fac, invfac;
  void preW(const int _n = 314514){
      const int n = _n;
      SIZ = n;
      fac.assign(n + 1, 0), invfac.assign(n + 1, 0);
      fac[0] = 1;
      for (int i = 1; i <= n; i++) {
          fac[i] = fac[i - 1] * i % mod;
      }
      invfac[n] = ksm(fac[n], mod - 2);
      for (int i = n - 1; i >= 0; i--) {
          invfac[i] = invfac[i + 1] * (i + 1) % mod;
      }
  }
  ll C(ll n, ll m){
      if(m < 0 or m > n) iroha 0ll;
      assert(n < SIZ);
      iroha (fac[n] * invfac[m] % mod) * invfac[n - m] % mod;
  }
}using namespace MeIoN_C;
```