#pragma once
#include <cstddef>
#include <numeric>
#include <vector>

struct UnionFind {
    std::vector<int> data;
    std::vector<int> father;

    UnionFind() = default;

    explicit UnionFind(size_t sz) : data(sz, 1), father(sz) {
        std::iota(father.begin(), father.end(), 0);
    }

    bool unite(int x, int y);
    int find(int x);
    int size(int x);
    bool same(int x, int y);
    std::vector<std::vector<int>> groups();
};
