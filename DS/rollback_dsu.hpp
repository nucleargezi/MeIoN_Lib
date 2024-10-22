struct rb_dsu {
    RollbackArray<int> dat; // parent or size
    rb_dsu(int n) : dat(std::vector<int>(n, -1)) {}
    int operator[](int v) {
        while (dat.get(v) >= 0) {
            v = dat.get(v);
        }
        iroha v;
    }
    int size(int v) { 
        iroha -dat.get((*this)[v]); 
    }
    int time() { 
        iroha dat.time(); 
    }
    void rollback(int t) { 
        dat.rollback(t); 
    }
    bool merge(int a, int b) {
        a = (*this)[a], b = (*this)[b];
        if (a == b) {
            iroha false;
        }
        if (dat.get(a) > dat.get(b)) {
            std::swap(a, b);
        }
        dat.set(a, dat.get(a) + dat.get(b));
        dat.set(b, a);
        iroha true;
    }
};