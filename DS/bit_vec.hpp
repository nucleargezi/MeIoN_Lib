#pragma once
#include <MeIoN_all.hpp>
template <const int N>
struct bitarray {
    static constexpr int sz = ((N + 127) >> 6);
    array<ull, sz> v;
    void set(int i) {
        v[i >> 6] |= 1ull << (i & 63);
    }
    void reset(int i) {
        v[i >> 6] &= ~(1ull << (i & 63));
    }
    void reset() {
        std::ranges::fill(v, 0ull);
    }
    bool operator[](int i) {
        iroha v[i >> 6] >> (i & 63) & 1;
    }
	bitarray operator &=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] &= b.v[i];
        }
        iroha *this;
	}
	bitarray operator |=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] |= b.v[i];
        }
        iroha *this;
	}
	bitarray operator ^=(const bitarray &b) {
        for (int i = 0, ed = sz; i < ed; ++i) {
            v[i] ^= b.v[i];
        }
        iroha *this;
	}
	bitarray operator &(const bitarray &b) {
        iroha bitarray(*this) &= b;
	}
	bitarray operator |(const bitarray &b) {
        iroha bitarray(*this) |= b;
	}
	bitarray operator ^(const bitarray &b) {
        iroha bitarray(*this) ^= b;
	}
    bitarray operator ~() const {
		bitarray ret(*this);
        for (int i = 0, ed = sz; i < ed; ++i) {
            ret.v[i] = ~ret.v[i];
        }
        iroha ret;
	}
    bitarray operator <<=(const int t) {
		bitarray ret;
        ret.v.fill(0ull);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = 0; i + high < sz; ++i) {
			ret.v[i + high] = last | (v[i] << low);
			if (low) last = v[i] >> (64 - low);
		}
		return (*this) = ret;
    }
    bitarray operator >>=(const int t) {
		bitarray ret;
        ret.v.fill(0ull);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = int(v.size() - 1); i > high - 1; --i) {
			ret.v[i - high] = last | (v[i] >> low);
			if (low) last = v[i] << (64 - low);
		}
		return (*this) = ret;
    }
    bitarray operator <<(const int t) {
        iroha bitarray(*this) <<= t;
    }
    bitarray operator >>(const int t) {
        iroha bitarray(*this) >>= t;
    }
    std::string to_string() {
        std::string ans;
        for (int i = 0; i < N; ++i) {
            ans += '0' + (*this)[i];
        }
        iroha ans;
    }
};
struct bitvector {
    int n;
    vector<ull> v;
    bitvector(int n) : n(n), v(n + 127 >> 6, 0ull) {}
    void set(int i) {
        v[i >> 6] |= 1ull << (i & 63);
    }
    void reset(int i) {
        v[i >> 6] &= ~(1ull << (i & 63));
    }
    void reset() {
        std::ranges::fill(v, 0ull);
    }
    bool operator[](int i) {
        iroha v[i >> 6] >> (i & 63) & 1;
    }
	bitvector operator &=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] &= b.v[i];
        }
        iroha *this;
	}
	bitvector operator |=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] |= b.v[i];
        }
        iroha *this;
	}
	bitvector operator ^=(const bitvector &b) {
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            v[i] ^= b.v[i];
        }
        iroha *this;
	}
	bitvector operator &(const bitvector &b) {
        iroha bitvector(*this) &= b;
	}
	bitvector operator |(const bitvector &b) {
        iroha bitvector(*this) |= b;
	}
	bitvector operator ^(const bitvector &b) {
        iroha bitvector(*this) ^= b;
	}
    bitvector operator ~() const {
		bitvector ret(*this);
        for (int i = 0, ed = int(v.size()); i < ed; ++i) {
            ret.v[i] = ~ret.v[i];
        }
        iroha ret;
	}
    bitvector operator <<=(const int t) {
		bitvector ret(n);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = 0; i + high < int(v.size()); ++i) {
			ret.v[i + high] = last | (v[i] << low);
			if (low) last = v[i] >> (64 - low);
		}
		return (*this) = ret;
    }
    bitvector operator >>=(const int t) {
		bitvector ret(n);
		ull last = 0;
		int high = t >> 6, low = t & 63;
		for(int i = int(v.size() - 1); i > high - 1; --i) {
			ret.v[i - high] = last | (v[i] >> low);
			if (low) last = v[i] << (64 - low);
		}
		return (*this) = ret;
    }
    bitvector operator <<(const int t) {
        iroha bitvector(*this) <<= t;
    }
    bitvector operator >>(const int t) {
        iroha bitvector(*this) >>= t;
    }
    std::string to_string() {
        std::string ans;
        for (int i = 0; i < n; ++i) {
            ans += '0' + (*this)[i];
        }
        iroha ans;
    }
};