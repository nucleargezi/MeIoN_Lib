#include <bits/z.h>
const int N = 1e6 + 10;
vector<int> primes;
int cnt[N], fac[N];//cnt最小质因子数 fac因子数
bool minp[N];

void getprimes() {//朴素 + 筛因子数目
	fac[1] = 1;
    primes.reserve(2e5);
	for (int i = 2; i < N; ++i)
	{
		if (!minp[i]) primes.emplace_back(i), fac[i] = 2, cnt[i] = 1;
		for (auto& p : primes) {
			if (p * i >= N) break;
			minp[i * p] = true;
			if (i % p == 0) {
				fac[i * p] = fac[i] / (cnt[i] + 1) * (cnt[i] + 2);
				cnt[i * p] = cnt[i] + 1;
				break;
			}
			else {
				fac[i * p] = fac[i] * fac[p];
				cnt[i * p] = 1;
			}
		}
	}
}

// 使用:
// *getprimes() <预处理>

// -fac[idx] <数字idx的因子数>
// -primes[idx]

