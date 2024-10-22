namespace MeIoN_namache{
    ll n;
    constexpr int working_sz = 1145141;
    int p[working_sz];
    string ss;
    void namache(string &s) {
        s = " " + s;
        ss = "";
        ss += '&';
        for (int i = 1; i <= n; i++) {
            ss += '#';
            ss += s[i];
        }
        ss += '#';
        ss += '*';
        for (int i = 1; i < ss.size(); i++) p[i] = 1;
        for (int i = 1, l = 1, r = 0; i + 1 < ss.size(); i++) {
            if (i <= r) p[i] = std::min(r - i + 1, p[l + r - i]);
            while (ss[i - p[i]] == ss[i + p[i]]) p[i] ++;
            if (i + p[i] - 1 > r) {
                l = i - p[i] + 1;
                r = i + p[i] - 1;
            }
        }
    }
    // [l, r]
    bool check(int l, int r) {
        int len = r - l + 1;
        l *= 2, r *= 2;
        return p[l + r >> 1] - 1 >= len;
    }
}using namespace MeIoN_namache;