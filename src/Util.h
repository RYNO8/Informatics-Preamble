#pragma once
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <deque>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <sstream>
#include "Constants.h"

namespace DS {
    /************************************************
     *                    DISPLAY                   *
     ************************************************/
    
    // Displays a pair:
    // (A, B)
    template<typename A, typename B> std::ostream& operator<<(std::ostream& out, const std::pair<A, B>& c) {
        out << '(' << c.first << ", " << c.second << ')';
        return out;
    }

    // Displays a vector:
    // [ T, T, ..., T ]
    template<typename T> std::ostream& operator<<(std::ostream& out, const std::vector<T>& v) {
        out << "[ ";
        for (int i = 0; i < (int)v.size(); ++i) {
            out << v[i];
            if (i != v.size() - 1) out << ", ";
        }
        out << " ]";
        return out;
    }

    // Displays an array:
    // T T T ... T
    #define printArr0(a, N) for (int i = 0; i < N; ++i) std::cout << a[i] << " ";
    #define printArr1(a, N) for (int i = 1; i <= N; ++i) std::cout << a[i] << " ";
    #define printArr(a, N) printArr0(a, N);

    // Displays a queue:
    // < T, T, ..., T >
    template<typename T> std::ostream& operator<<(std::ostream& out, const std::queue<T>& q) {
        out << "< ";
        for (int i = 0; i < (int)q.size(); ++i) {
            out << q.front();
            q.push(q.front());
            q.pop();
            if (i != q.size() - 1) out << ", ";
        }
        out << " >";
        return out;
    }

    // Displays a stack:
    // < T, T, ..., T >
    template<typename T> std::ostream& operator<<(std::ostream& out, const std::stack<T>& s) {
        out << "< ";
        for (int i = 0; i < (int)s.size(); ++i) {
            out << s.top();
            s.push(s.top());
            s.pop();
            if (i != s.size() - 1) out << ", ";
        }
        out << " >";
        return out;
    }

    // Displays a deque:
    // < T, T, ..., T >
    template<typename T> std::ostream& operator<<(std::ostream& out, const std::deque<T>& q) {
        out << "< ";
        for (int i = 0; i < (int)q.size(); ++i) {
            out << q.front();
            q.push_back(q.front());
            q.pop_front();
            if (i != q.size() - 1) out << ", ";
        }
        out << " >";
        return out;
    }

    // Displays a map:
    // { A:B, A:B, ..., A:B }
    template<typename A, typename B> std::ostream& operator<<(std::ostream& out, const std::map<A, B> &m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << it->first << ':' << it->second;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays an unordered map:
    // { A:B, A:B, ..., A:B }
    template<typename A, typename B> std::ostream& operator<<(std::ostream& out, const std::unordered_map<A, B>& m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << it->first << ':' << it->second;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays a multimap:
    // { A:B, A:B, ..., A:B }
    template<typename A, typename B> std::ostream& operator<<(std::ostream& out, const std::multimap<A, B> &m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << it->first << ':' << it->second;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays an unordered multimap:
    // { A:B, A:B, ..., A:B }
    template<typename A, typename B> std::ostream& operator<<(std::ostream& out, const std::unordered_multimap<A, B>& m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << it->first << ':' << it->second;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays a set:
    // { T, T, ..., T }
    template<typename T> std::ostream& operator<<(std::ostream& out, const std::set<T>& m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << *it;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays an unordered set:
    // { T, T, ..., T }
    template<typename T> std::ostream & operator<<(std::ostream & out, const std::unordered_set<T> & m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << *it;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays a multiset:
    // { T, T, ..., T }
    template<typename T> std::ostream & operator<<(std::ostream & out, const std::multiset<T> & m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << *it;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    // Displays an unordered multiset:
    // { T, T, ..., T }
    template<typename T> std::ostream & operator<<(std::ostream & out, const std::unordered_multiset<T> & m) {
        out << "{ ";
        for (auto it = m.begin(); it != m.end(); ) {
            out << *it;
            if (++it != m.end()) out << ", ";
        }
        out << " }";
        return out;
    }

    /************************************************
     *           NUMERICAL TYPE UTILITIES           *
     ************************************************/

    // O(log a + log b)
    // @note `a` and `b` should be non-negative?
    template<typename T> T gcd(T x, T y) {
        while (y != 0) {
            T temp = y;
            y = x % temp;
            x = temp;
        }
        return x;
    }

    // O(log a + log b)
    // @note May overflow T
    template<typename T> T lcm(T a, T b) {
        return (a / gcd(a, b)) * b;
    }

    // O(12)
    // @returns the number of set bits in the binary represetnation of `x`
    uint popcount(uint x) {
        // note: faster than builtin __builtin_popcount
        x -= (x >> 1) & 0x55555555;
        x = (x & 0x33333333) + ((x >> 2) & 0x333333333);
        x = (x + (x >> 4)) & 0xf0f0f0f;
        return (x * 0x1010101) >> 24;
    }
    int popcount(int x) {
        assert(x >= 0);
        return popcount((uint)x);
    }
    ull popcount(ull x) {
        // note: faster than builtin __builtin_popcountll
        x -= (x >> 1) & 0x5555555555555555;
        x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333);
        x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;
        return (x * 0x0101010101010101) >> 56;
    }
    ll popcount(ll x) {
        assert(x >= 0);
        return popcount((ull)x);
    }

    // O(1)
    // @returns Whether the number is negative (-1), zero (0) or positive (+1)
    template<typename T> T sgn(T x) {
        return T(x > 0) - T(x < 0);
    }

    // O(1)
    // Sets `a` to the minimum of `a` and `b`
    // @note `a` is provided by reference
    template<typename T> void pMin(T& a, T b) {
        if (b < a) a = b;
    }

    // O(1)
    // Sets `a` to the maximum of `a` and `b`
    // @note `a` is provided by reference
    template<typename T> void pMax(T& a, T b) {
        if (b > a) a = b;
    }

    // O(N)
    // @returns The number of characters of the printed representation of `obj`
    template<typename T> int reprLen(T obj) {
        std::stringstream s;
        s << obj;
        return s.str().size();
    }

    /************************************************
     *                OTHER UTILITIES               *
     ************************************************/

    double timeNow() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
        ).count() / 1000.;
    }
};