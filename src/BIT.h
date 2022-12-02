#pragma once
#include <assert.h>
#include <iostream>
#include <vector>
#include "Constants.h"

namespace DS {
    // TODO: initialisation from vector, array
    // TODO: function to fill range with value
    // TODO: conversion between Grid

    enum BIT_mode {
        RANGE_QUERY_POINT_UPDATE = 0,
        POINT_QUERY_RANGE_UPDATE = 1
    };

    // T() is left and right identity
    // T operator+(T, T) is associative
    // T operator-(T, T) is inverse
    template <class T, BIT_mode mode, int ...Ns>
    class BIT {
private:
        T val = T();

public:
        size_t size() const {
            return 1;
        }

        size_t dimensions() const {
            return 0;
        }

        std::vector<size_t> shape() const {
            return {};
        }

        void addIndex(T v) {
            val += v;
        }

        void subtractIndex(T v) {
            val -= v;
        }

        void addRange(T v) {
            val += v;
        }

        void subtractRange(T v) {
            val -= v;
        }

        T queryIndex() const {
            return val;
        }
        
        T querySum() const {
            return val;
        }
    };

    template<class T, BIT_mode mode, int N, int... Ns>
    class BIT<T, mode, N, Ns...> {
private:
        BIT<T, mode, Ns...> bit[N + 1];

public:
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

        
        template<typename... Args> void addIndex(int pos, Args... args) {
            assert(0 <= pos && pos < N && "index out of range");
            if (mode == RANGE_QUERY_POINT_UPDATE) {
                for (++pos; pos <= N; pos += pos & -pos) bit[pos].addIndex(args...);
            } else if (mode == POINT_QUERY_RANGE_UPDATE) {
                for (int x = pos + 1; x; x -= x & -x) bit[x].addIndex(args...);
                for (int x = pos + 0; x; x -= x & -x) bit[x].subtractIndex(args...);
            }
        }
        
        template<typename... Args> void subtractIndex(int pos, Args... args) {
            assert(0 <= pos && pos < N && "index out of range");
            assert(mode == POINT_QUERY_RANGE_UPDATE);
            for (int x = pos + 1; x; x -= x & -x) bit[x].subtractIndex(args...);
            for (int x = pos + 0; x; x -= x & -x) bit[x].addIndex(args...);
        }

        template<typename... Args> void addRange(int l, int r, Args... args) {
            assert(0 <= l && l <= r && r < N && "index out of range");
            assert(mode == POINT_QUERY_RANGE_UPDATE);
            for (int x = r + 1; x; x -= x & -x) bit[x].addRange(args...);
            for (int x = l + 0; x; x -= x & -x) bit[x].subtractRange(args...);
        }
        
        template<typename... Args> void subtractRange(int l, int r, Args... args) {
            assert(0 <= l && l <= r && r < N && "index out of range");
            assert(mode == POINT_QUERY_RANGE_UPDATE);
            for (int x = r + 1; x; x -= x & -x) bit[x].subtractRange(args...);
            for (int x = l + 0; x; x -= x & -x) bit[x].addRange(args...);
        }

        template<typename... Args> T queryIndex(int pos, Args... args) const {
            assert(0 <= pos && pos < N && "index out of range");
            if (mode == RANGE_QUERY_POINT_UPDATE) {
                T res = 0;
                for (int x = pos + 1; x; x -= x & -x) res += bit[x].queryIndex(args...);
                for (int x = pos + 0; x; x -= x & -x) res -= bit[x].queryIndex(args...);
                return res;
            } else if (mode == POINT_QUERY_RANGE_UPDATE) {
                T res = 0;
                for (++pos; pos <= N; pos += pos & -pos) res += bit[pos].queryIndex(args...);
                return res;
            }
        }
        
        template<typename... Args> T querySum(int l, int r, Args... args) const {
            assert(0 <= l && l <= r && r < N && "index out of range");
            assert(mode == RANGE_QUERY_POINT_UPDATE);
            T res = 0;
            for (int x = r + 1; x; x -= x & -x) res += bit[x].querySum(args...);
            for (int x = l + 0; x; x -= x & -x) res -= bit[x].querySum(args...);
            return res;
        }
    };

    template<class T, BIT_mode mode> std::ostream& operator<<(std::ostream& out, const BIT<T, mode> b) {
        out << "[ " << b.queryIndex() << " ]\n";
        return out;
    }

    template<class T, BIT_mode mode, int N> std::ostream& operator<<(std::ostream& out, const BIT<T, mode, N> b) {
        out << "[ ";
        for (int i = 0; i < N; ++i) out << b.queryIndex(i) << ' ';
        out << " ]\n";
        return out;
    }

    template<class T, BIT_mode mode, int N, int M> std::ostream& operator<<(std::ostream& out, const BIT<T, mode, N, M> b) {
        out << "[ \n";
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < M; ++j) out << ' ' << b.queryIndex(i, j);
            out << '\n';
        }
        out << "]\n";
        return out;
    }
};