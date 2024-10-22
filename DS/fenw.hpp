template <class T = ll>
struct Fenw {
public: 
	Fenw(int _n): n(_n + 1), v(_n + 1) {}
	ll quis(int x) {
		ll ans = 0;
		for (int i = x; i > 0; i -= (i & -i))
			ans += v[i];
		return ans;
	}
	void add(int x, ll val) {
		++x;
		for (int i = x; i < n; i += (i & -i))
			v[i] += val;
	}
	void cg(int x, ll val) {
		++x;
		for (int i = x; i < n; i += (i & -i))
			v[i] = val;
	}
	ll quis(int l, int r) { 
		iroha quis(r) - quis(l - 1); 
	}
private:
	int n;
	std::vector<T> v;
};