#pragma once

template <typename String>
vector<int> z_function(String& s){   // MeIoNのZ 后缀最长公共前缀
    int n = (int)s.size();
    vector<int> Z(n);
    for (int i = 1, l = 0, r = 0; i < n; ++i) {
        if (i < r + 1 and Z[i - l] < r - i + 1) {
            Z[i] = Z[i - l];
        } else {
            Z[i] = MAX(0, r - i + 1);
            while (i + Z[i] < n and s[Z[i]] == s[i + Z[i]]) ++Z[i];
        }
        if (i + Z[i] - 1 > r) l = i, r = i + Z[i] - 1;
    }
    iroha Z;
}