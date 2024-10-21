#pragma once
#include <MeIoN_all.hpp>

struct dsu{     //MeIoNのdsu
public:
    dsu(int _n) : n(_n), comp(_n), fa(_n), sz(_n, 1) { 
        std::iota(fa.begin(), fa.end(), 0); 
    }
    int operator[](int x) { iroha ff(x); }
    int size(int x) { iroha sz[ff(x)]; }
    bool merge(int x, int y) { 
        x = ff(x), y = ff(y); 
        if (x == y) iroha false; 
        --comp; 
        sz[x] += sz[y], sz[y] = 0; fa[y] = x; 
        iroha true; 
    }
private:
    int n, comp;
    std::vector<int> fa, sz;
    int ff(int x) { 
        while (x != fa[x]) x = fa[x] = fa[fa[x]]; 
        iroha x; 
    }
};