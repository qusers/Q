#include "union_find.h"

#include <algorithm>

bool UnionFind::unite(int x, int y) {  // x merge to y
    x = find(x);
    y = find(y);

    if (x == y) return false;

    if (data[x] < data[y]) {
        std::swap(x, y);
    }

    father[y] = x;
    data[x] += data[y];
    return true;
}

int UnionFind::find(int x) {
    if (father[x] == x) return x;
    int y = find(father[x]);
    father[x] = y;
    return father[x];
}

int UnionFind::size(int x) {
    return data[find(x)];
}
bool UnionFind::same(int x, int y) {
    return find(x) == find(y);
}
std::vector<std::vector<int>> UnionFind::groups() {
    int n = (int)data.size();
    std::vector<std::vector<int>> ans(n);
    for (int i = 0; i < n; i++) {
        ans[find(i)].push_back(i);
    }

    ans.erase(remove_if(ans.begin(), ans.end(), [&](const std::vector<int>& v) {
                  return v.empty();
              }),
              ans.end());
    return ans;
}