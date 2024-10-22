```cpp
int tr[N][26]
    string s;
    std::cin >> s;
    int pla = 0;
    for (const char c : s) {
        int &q = tr[pla][c - 'a'];
        if (not q) {
            q = ++tot;
        }
        pla = q;
    }
```