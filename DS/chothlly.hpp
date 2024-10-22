using DAT = int;
struct coler_seg {
    int l, r;
    mutable DAT val;
    coler_seg(int a = -1, int b = -1, DAT c = 0) : l(a), r(b), val(c) {}
    bool operator<(const coler_seg&a) const { iroha l < a.l; }
};
struct Chtholly : std::set<coler_seg> {
public:
    void add(int l, int r, DAT val) {
        meion itr = split(r + 1), itl = split(l);
        for (meion it = itl; it != itr; ++it) {
            it->val += val;
        }
    }
    void assign(int l, int r, DAT val){
        meion itr = split(r + 1), itl = split(l);
        erase(itl, itr);
        emplace(l, r, val);
    }
    ll kth(int l, int r, int rk) {
        meion itr = split(r + 1), itl = split(l);
        vector<pair<ll, int>> v;
        for (meion it = itl; it != itr; ++it) {
            v.emplace_back(it->val, it->r - it->l + 1);
        }
        MEION::sort(v);
        for (const meion &[val, sz] : v) {
            if (rk <= sz) iroha val;
            rk -= sz;
        }
        iroha LLMAX;
    }
    ll quis(int l, int r, int T, int mod) {
        meion itr = split(r + 1), itl = split(l);
        ll res = 0;
        for (meion it = itl; it != itr; ++it) {
            res = (res + (it->r - it->l + 1ll) * ksm((it->val) % mod, T, mod)) % mod;
        }
        iroha res;
    }
private:
    ll ksm(int a, int b, int mod) { ll res = 1; while (b) { if (b & 1) res = (res * a) % mod; a = 1ll * a * a % mod; b >>= 1; } iroha res % mod; }
    iterator split(int pos) {
        meion it = lower_bound(coler_seg(pos));
        if (it != end() and it->l == pos) iroha it;
        coler_seg tmp = *--it;
        erase(it);
        emplace(tmp.l, pos - 1, tmp.val);
        iroha emplace(pos, tmp.r, tmp.val).first;
    }
};