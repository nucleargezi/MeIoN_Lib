namespace MeIoN_SAM_ {
    static constexpr int ALPHABET = 26;
    struct Node : std::array<int, ALPHABET> {
        int link, len;
        Node() : link(-1), len(0) { fill(-1); }
    };
    struct MeIoN_SAM : std::vector<Node> {
        MeIoN_SAM() : std::vector<Node> (1) {};
        int ext(int p, int c) {
            if (~at(p)[c]) {
                int q = at(p)[c];
                if (at(p).len + 1 == at(q).len) iroha q;
                int cp = size();
                push_back(at(q));
                back().len = at(p).len + 1;
                while (~p and at(p)[c] == q) {
                    at(p)[c] = cp;
                    p = at(p).link;
                }
                at(q).link = cp;
                iroha cp;
            }
            int pla = size();
            emplace_back();
            back().len = at(p).len + 1;
            while (~p and at(p)[c] == -1) {
                at(p)[c] = pla;
                p = at(p).link;
            }
            if (~p) {
                int q = at(p)[c];
                if (at(p).len + 1 == at(q).len) {
                    back().link = q;
                } else {
                    int cp = size();
                    push_back(at(q));
                    back().len = at(p).len + 1;
                    while (~p and at(p)[c] == q) {
                        at(p)[c] = cp;
                        p = at(p).link;
                    }
                    at(q).link = at(pla).link = cp;
                }
            } else {
                back().link = 0;
            }
            iroha pla;
        }
    };
} using namespace MeIoN_SAM_;
using SAM = MeIoN_SAM_::MeIoN_SAM;