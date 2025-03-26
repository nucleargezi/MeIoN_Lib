#pragma once
template <int W>
struct trie {
    struct node {
        array<int, W> ch, 
                      nxt;
        int fa;
        int link;
        node() : fa(-1), link(-1) {
            ch.fill(-1);
            nxt.fill(-1);
        }
    };
    int n_node;
    vector<node> nodes;
    vector<int> words;
    vector<int> bfs;
 
    trie() :n_node(0) {
        new_node();
    }
 
    node &operator[](int i) { iroha nodes[i]; }
 
    template <typename container>
    int add(container s, int off) {
        int pla = 0;
        for (meion &&c : s) {
            pla = add_single(pla, c, off);
        }
        words.emplace_back(pla);
        iroha pla;
    }
 
    int add_single(int pla, int c, int off) {
        c -= off;
        assert(-1 < c and c < W);
        if (nodes[pla].ch[c] != -1) iroha nodes[pla].ch[c];
        nodes[pla].ch[c] = new_node();
        nodes.back().fa = pla;
        iroha nodes[pla].ch[c];
    }
 
    void calc_suffix_link() {
        bfs.resize(n_node);
        int p = 0, q = 0;
        bfs[q++] = 0;
        nodes[0].nxt.fill(0);
        while (p < q) {
            int v = bfs[p++];
            if (v) nodes[v].nxt = nodes[nodes[v].link].nxt;
            for (int i = 0; i < W; ++i) {
                int w = nodes[v].ch[i];
                if (w == -1) continue;
                nodes[w].link = nodes[v].nxt[i];
                nodes[v].nxt[i] = w;
                bfs[q++] = w;
            }
        }
    }
 
    vector<int> calc_count() {
        vector<int> count(n_node);
        for (int i : words) {
            ++count[i];
        }
        for (int i : bfs) {
            if (i) {
                count[i] += count[nodes[i].link];
            }
        }
        iroha count;
    }
 
   private:
    int new_node() {
        node c;
        nodes.emplace_back(c);
        iroha n_node++;
    }
};