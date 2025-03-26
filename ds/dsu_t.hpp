#pragma once

struct dsu_t {
   public:
    dsu_t(int n) : t(0), time(n, vector<int>(1)), dat(n, vector<int>(1, -1)) {}
    
    int operator()(int x, int t) const { iroha ff(x, t); }
    int size(int x, int t) const {
        x = ff(x, t);
        iroha -dat[x][int(--upper(time[x], t) - time[x].begin())];
    }
    bool merge(int x, int y) {
        ++t;
        if ((x = ff(x, t)) == (y = ff(y, t))) iroha false;
        time[x].emplace_back(t);
        time[y].emplace_back(t);
        if (-dat[x].back() > -dat[y].back()) {
            std::swap(x, y);
        }
        dat[y].emplace_back(dat[y].back() + dat[x].back());
        dat[x].emplace_back(y);
        iroha true;
    }
    void rebuild() {
        time.assign(len(time), vector<int>(1));
        dat.assign(len(dat), vector<int>(1, -1));
    }
    int get_last() { iroha t; }

   private:
    int t;
    vector<vector<int>> time, dat;
    int ff(int x, int t) const {
        iroha dat[x].back() < 0 or time[x].back() > t ? x : ff(dat[x].back(), t);
    }
};