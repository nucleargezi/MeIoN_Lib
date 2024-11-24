#pragma once
template <bool dir = false>
void read_graph(vector<vector<int>> &v, const int edges) {
    if constexpr (dir) {
        for (int i = 0, x, y; i < edges; ++i) {
            std::cin >> x >> y, --x, --y;
            v[x].emplace_back(y);
        }
    } else {
        for (int i = 0, x, y; i < edges; ++i) {
            std::cin >> x >> y, --x, --y;
            v[x].emplace_back(y);
            v[y].emplace_back(x);
        }
    }
}
template <bool dir = false, typename T>
void read_graph_w(vector<vector<pair<int, T>>> &v, const int edges) {
    if constexpr (dir) {
        for (int i = 0, x, y; i < edges; ++i) {
            T w;
            std::cin >> x >> y >> w;
            v[x].emplace_back(y, w);
        }
    } else {
        for (int i = 0, x, y; i < edges; ++i) {
            T w;
            std::cin >> x >> y >> w;
            v[x].emplace_back(y, w);
            v[y].emplace_back(x, w);
        }
    }
}
template <bool dir = false>
void read_tree(vector<vector<int>> &v) {
    const int edges = v.size();
    if constexpr (dir) {
        for (int i = 1, x, y; i < edges; ++i) {
            std::cin >> x >> y, --x, --y;
            v[x].emplace_back(y);
        }
    } else {
        for (int i = 1, x, y; i < edges; ++i) {
            std::cin >> x >> y, --x, --y;
            v[x].emplace_back(y);
            v[y].emplace_back(x);
        }
    }
}
template <bool dir = false, typename T>
void read_tree_w(vector<vector<pair<int, T>>> &v) {
    const int n = v.size();
    if constexpr (dir) {
        for (int i = 1, x, y; i < n; ++i) {
            T w;
            std::cin >> x >> y >> w;
            v[x].emplace_back(y, w);
        }
    } else {
        for (int i = 1, x, y; i < n; ++i) {
            T w;
            std::cin >> x >> y >> w;
            v[x].emplace_back(y, w);
            v[y].emplace_back(x, w);
        }
    }
}