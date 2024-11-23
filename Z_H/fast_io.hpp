#pragma once
namespace fast_io {
    static constexpr uint32_t SZ = 1 << 17;
    char ibuf[SZ];
    char obuf[SZ];
    char out[100];
    // pointer of ibuf, obuf
    uint32_t pil = 0, pir = 0, por = 0;

    struct Pre {
        char num[10000][4];
        constexpr Pre() : num() {
            for (int i = 0; i < 10000; i++) {
                int n = i;
                for (int j = 3; j >= 0; j--) {
                    num[i][j] = n % 10 | '0';
                    n /= 10;
                }
            }
        }
    } constexpr pre;

    inline void load() {
        memcpy(ibuf, ibuf + pil, pir - pil);
        pir = pir - pil + fread(ibuf + pir - pil, 1, SZ - pir + pil, stdin);
        pil = 0;
        if (pir < SZ) ibuf[pir++] = '\n';
    }

    inline void flush() {
        fwrite(obuf, 1, por, stdout);
        por = 0;
    }
    void rd(char &c) {
        do {
            if (pil + 1 > pir) load();
            c = ibuf[pil++];
        } while (isspace(c));
    }

    void rd(string &x) {
        x.clear();
        char c;
        do {
            if (pil + 1 > pir) load();
            c = ibuf[pil++];
        } while (isspace(c));
        do {
            x += c;
            if (pil == pir) load();
            c = ibuf[pil++];
        } while (!isspace(c));
    }

    template <typename T>
    void rd_real(T &x) {
        string s;
        rd(s);
        x = stod(s);
    }

    template <typename T>
    void rd_integer(T &x) {
        if (pil + 100 > pir) load();
        char c;
        do c = ibuf[pil++];
        while (c < '-');
        bool minus = 0;
        if constexpr (std::is_signed<T>::value || std::is_same_v<T, i128>) {
            if (c == '-') {
                minus = 1, c = ibuf[pil++];
            }
        }
        x = 0;
        while ('0' <= c) {
            x = x * 10 + (c & 15), c = ibuf[pil++];
        }
        if constexpr (std::is_signed<T>::value || std::is_same_v<T, i128>) {
            if (minus) x = -x;
        }
    }

    void rd(int &x) { rd_integer(x); }
    void rd(ll &x) { rd_integer(x); }
    void rd(i128 &x) { rd_integer(x); }
    void rd(uint &x) { rd_integer(x); }
    void rd(ull &x) { rd_integer(x); }
    void rd(u128 &x) { rd_integer(x); }
    void rd(double &x) { rd_real(x); }
    void rd(long double &x) { rd_real(x); }
    void rd(f128 &x) { rd_real(x); }

    template <class T, class U>
    void rd(pair<T, U> &p) {
        return rd(p.first), rd(p.second);
    }
    template <size_t N = 0, typename T>
    void rd_tuple(T &t) {
        if constexpr (N < std::tuple_size<T>::value) {
            auto &x = std::get<N>(t);
            rd(x);
            rd_tuple<N + 1>(t);
        }
    }
    template <class... T>
    void rd(std::tuple<T...> &tpl) {
        rd_tuple(tpl);
    }

    template <size_t N = 0, typename T>
    void rd(array<T, N> &x) {
        for (auto &d : x) rd(d);
    }
    template <class T>
    void rd(vector<T> &x) {
        for (auto &d : x) rd(d);
    }

    void read() {}
    template <class H, class... T>
    void read(H &h, T &...t) {
        rd(h), read(t...);
    }

    void wt(const char c) {
        if (por == SZ) flush();
        obuf[por++] = c;
    }
    void wt(const string s) {
        for (char c : s) wt(c);
    }
    void wt(const char *s) {
        size_t len = strlen(s);
        for (size_t i = 0; i < len; i++) wt(s[i]);
    }

    template <typename T>
    void wt_integer(T x) {
        if (por > SZ - 100) flush();
        if (x < 0) {
            obuf[por++] = '-', x = -x;
        }
        int outi;
        for (outi = 96; x >= 10000; outi -= 4) {
            memcpy(out + outi, pre.num[x % 10000], 4);
            x /= 10000;
        }
        if (x >= 1000) {
            memcpy(obuf + por, pre.num[x], 4);
            por += 4;
        } else if (x >= 100) {
            memcpy(obuf + por, pre.num[x] + 1, 3);
            por += 3;
        } else if (x >= 10) {
            int q = (x * 103) >> 10;
            obuf[por] = q | '0';
            obuf[por + 1] = (x - q * 10) | '0';
            por += 2;
        } else
            obuf[por++] = x | '0';
        memcpy(obuf + por, out + outi + 4, 96 - outi);
        por += 96 - outi;
    }

    template <typename T>
    void wt_real(T x) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(15) << double(x);
        std::string s = oss.str();
        wt(s);
    }

    void wt(int x) { wt_integer(x); }
    void wt(ll x) { wt_integer(x); }
    void wt(i128 x) { wt_integer(x); }
    void wt(uint x) { wt_integer(x); }
    void wt(ull x) { wt_integer(x); }
    void wt(u128 x) { wt_integer(x); }
    void wt(double x) { wt_real(x); }
    void wt(long double x) { wt_real(x); }
    void wt(f128 x) { wt_real(x); }
    void wt(std::ios_base &(*__pf)(std::ios_base &)) {}
    void wt(const std::_Setprecision &x) {}

    template <class T, class U>
    void wt(const pair<T, U> val) {
        wt(val.first);
        wt(' ');
        wt(val.second);
    }
    template <size_t N = 0, typename T>
    void wt_tuple(const T t) {
        if constexpr (N < std::tuple_size<T>::value) {
            if constexpr (N > 0) {
                wt(' ');
            }
            const auto x = std::get<N>(t);
            wt(x);
            wt_tuple<N + 1>(t);
        }
    }
    template <class... T>
    void wt(std::tuple<T...> tpl) {
        wt_tuple(tpl);
    }
    template <class T, size_t S>
    void wt(const array<T, S> val) {
        auto n = val.size();
        for (size_t i = 0; i < n; i++) {
            if (i) wt(' ');
            wt(val[i]);
        }
    }
    template <class T>
    void wt(const vector<T> val) {
        auto n = val.size();
        for (size_t i = 0; i < n; i++) {
            if (i) wt(' ');
            wt(val[i]);
        }
    }
    // gcc expansion. called automaticall after main.
    void __attribute__((destructor)) _d() { flush(); }
    
    struct io_auxiliary {
        void sync_with_stdio(bool ok) { }
    } *auxiliary_io;

    struct meion_fast_io {
        template<typename T>
        friend meion_fast_io& operator>>(meion_fast_io &in, T &c) {
            rd(c);
            iroha in;
        }
        template<typename T>
        friend meion_fast_io& operator<<(meion_fast_io &out, const T &c) {
            wt(c);
            iroha out;
        }
        io_auxiliary *tie(std::nullptr_t x) {
            iroha auxiliary_io;
        }
    } fin, fout;
}
#define fast
#define cin fin
#define cout fout
namespace std {
    using fast_io::cin, fast_io::cout;
}