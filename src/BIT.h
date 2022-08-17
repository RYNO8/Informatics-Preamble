#pragma once
#include "Constants.h"

// TODO: initialisation from vector, array
// TODO: function to fill range with value
// TODO: conversion between Grid

template <class T, int ...Ns> struct BIT {
    T val = T(0);

    size_t size() const {
        return 1;
    }

    size_t dimensions() const {
        return 0;
    }

    void add(T v) {
        val += v;
    }

    T getSum() const {
        return val;
    }

    T getVal() const {
        return val;
    }
};

template<class T, int N, int... Ns> struct BIT<T, N, Ns...> {
    BIT<T, Ns...> bit[N + 1];

    size_t size() const {
        return N * bit[0].size();
    }

    size_t dimensions() const {
        return 1 + bit[0].dimensions();
    }

    std::vector<size_t> shape() const {
        std::vector<size_t> output = { N };
        std::vector<size_t> inner = bit[0].shape();
        output.insert(output.end(), inner.begin(), inner.end());
        return output;
    }

    template<typename... Args> void add(int pos, Args... args) {
        assert(0 <= pos && pos < N && "index out of range");
        for (++pos; pos <= N; pos += pos & -pos) bit[pos].add(args...);
    }

    template<typename... Args> T getSum(int l, int r, Args... args) const {
        assert(0 <= l && l <= r && r < N && "index out of range");
        T res = 0;
        for (int x = r + 1; x; x -= x & -x) res += bit[x].getSum(args...);
        for (int x = l + 0; x; x -= x & -x) res -= bit[x].getSum(args...);
        return res;
    }

    template<typename... Args> T getVal(int pos, Args... args) const {
        assert(0 <= pos && pos < N && "index out of range");
        T res = 0;
        for (int x = pos + 1; x; x -= x & -x) res += bit[x].getVal(args...);
        for (int x = pos + 0; x; x -= x & -x) res -= bit[x].getVal(args...);
        return res;
    }
};

template<class T> std::ostream& operator<<(std::ostream& out, const BIT<T> b) {
    out << "[ " << b.val << " ]\n";
    return out;
}

template<class T, int N> std::ostream& operator<<(std::ostream& out, const BIT<T, N> b) {
    out << "[ ";
    for (int i = 0; i < N; ++i) out << b.getVal(i) << ' ';
    out << " ]\n";
    return out;
}

template<class T, int N, int M> std::ostream& operator<<(std::ostream& out, const BIT<T, N, M> b) {
    out << "[ \n";
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) out << ' ' << b.getVal(i, j);
        out << '\n';
    }
    out << "]\n";
    return out;
}