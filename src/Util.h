#ifndef UTIL_H
#define UTIL_H
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <utility>
#include <vector>
#include <queue>
#include <stack>
#include <deque>
#include <map>
#include <set>
#include "Constants.h"


namespace std {
    template<>
    struct hash<pair<int, int>> {
        size_t operator()(const pair<int, int> k) const {
            return hash<long long>()(((long long)k.first) << 32 | (long long) k.second);
        }
    };

    template<>
    struct hash<pair<unsigned int, unsigned int>> {
        size_t operator()(const pair<unsigned int, unsigned int> k) const {
            return hash<unsigned long long>()(((unsigned long long)k.first) << 32 | (unsigned long long) k.second);
        }
    };
}

namespace DS {
    /************************************************
     *                    DISPLAY                   *
     ************************************************/
    
    // Displays a pair:
    // (A, B)
    template<typename A, typename B>
    std::ostream& operator<<(std::ostream &out, const std::pair<A, B> &c) {
        out << '(' << c.first << ", " << c.second << ')';
        return out;
    }

    // Displays an array:
    // [ T T T ... T ]
    #define printArr0(a, N) std::cout << "[ "; for (int i = 0; i < N; ++i, std::cout << (i < N ? ", " : "")) std::cout << a[i]; std::cout << " ]";
    #define printArr1(a, N) std::cout << "[ "; for (int i = 1; i <= N; ++i, std::cout << (i < N ? ", " : "")) std::cout << a[i]; std::cout << " ]";
    #define printArr(a, N) printArr0(a, N)

    template<typename T> struct capture {
        T begin, end;
    };
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const capture<T> &range) {
        out << "[ ";
        for (T it = range.begin; it != range.end; ) {
            out << *it;
            if (++it != range.end) out << " ";
        }
        out << " ]";
        return out;
    }

    // Displays a vector:
    // [ T T ... T ]
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const std::vector<T> &v) {
        out << "[ ";
        // copy(v.cbegin(), v.cend(), std::ostream_iterator<T>(out, " "));
        for (auto val : v) out << val << ' ';
        out << ']';
        return out;
    }

    // Displays a std::array
    // [ T, T, ..., T ]
    template<typename T, std::size_t N>
    std::ostream& operator<<(std::ostream& out, const std::array<T, N> &arr){
        out << "[ ";
        // copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(out, " "));
        for (auto val : arr) out << val << ' ';
        out << ']';
        return out;
    }

    // Displays a queue:
    // < T, T, ..., T > (In the order they would be popped)
    // @note takes copy of queue
    template<typename T>
    std::ostream& operator<<(std::ostream &out, std::queue<T> q) {
        out << "<";
        while (!q.empty()) {
            out << ' ' << q.front();
            q.pop();
            if (!q.empty()) out << ',';
        }
        out << " >";
        return out;
    }

    // Displays a stack:
    // < T, T, ..., T > (In the order they would be popped)
    // @note takes copy of stack
    template<typename T>
    std::ostream& operator<<(std::ostream &out, std::stack<T> s) {
        out << '<';
        while (!s.empty()) {
            out << ' ' << s.top();
            s.pop();
            if (!s.empty()) out << ',';
        }
        out << " >";
        return out;
    }

    // Displays a deque:
    // < T, T, ..., T >
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const std::deque<T>& dq) {
        out << '<';
        for (auto it = dq.begin(); it != dq.end(); ) {
            out << ' ' << *it;
            if (++it != dq.end()) out << ',';
        }
        out << " >";
        return out;
    }

    // Displays a map:
    // { A: B, A: B, ..., A: B }
    template<typename A, typename B>
    std::ostream& operator<<(std::ostream &out, const std::map<A, B> &m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << it->first << ": " << it->second;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays an unordered map:
    // { A: B, A: B, ..., A: B }
    template<typename A, typename B>
    std::ostream& operator<<(std::ostream &out, const std::unordered_map<A, B>& m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << it->first << ": " << it->second;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays a multimap:
    // { A: B, A: B, ..., A: B }
    template<typename A, typename B>
    std::ostream& operator<<(std::ostream &out, const std::multimap<A, B> &m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << it->first << ": " << it->second;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays an unordered multimap:
    // { A: B, A: B, ..., A: B }
    template<typename A, typename B>
    std::ostream& operator<<(std::ostream &out, const std::unordered_multimap<A, B>& m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << it->first << ": " << it->second;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays a set:
    // { T, T, ..., T }
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const std::set<T>& m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << *it;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays an unordered set:
    // { T, T, ..., T }
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const std::unordered_set<T> &m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << *it;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays a multiset:
    // { T, T, ..., T }
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const std::multiset<T> &m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << *it;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }

    // Displays an unordered multiset:
    // { T, T, ..., T }
    template<typename T>
    std::ostream& operator<<(std::ostream &out, const std::unordered_multiset<T> &m) {
        out << '{';
        for (auto it = m.begin(); it != m.end(); ) {
            out << ' ' << *it;
            if (++it != m.end()) out << ',';
        }
        out << " }";
        return out;
    }


    /************************************************
     *           NUMERICAL TYPE UTILITIES           *
     ************************************************/

    // O(log a + log b)
    // @note `a` and `b` should be non-negative?
    template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
    T gcd(T x, T y) {
        while (y != 0) {
            T temp = y;
            y = x % temp;
            x = temp;
        }
        return x;
    }

    // O(log a + log b)
    // @note May overflow T
    template<typename T, std::enable_if_t<std::is_integral<T>::value, bool> = true>
    T lcm(T a, T b) {
        return (a / gcd(a, b)) * b;
    }

    // O(log N)
    // extended euclidean algorithm
    // @returns gcd(a, b)
    // sets x and y such that a*x + b*y = gcd(a, b)
    template<
        typename U, std::enable_if_t<std::is_integral<U>::value, bool> = true,
        typename S, std::enable_if_t<std::is_integral<S>::value, bool> = true
    >
    U extendedEuclidean(U a, U b, S &x, S &y) {
        if (b == 0) {
            x = S(1);
            y = S(0);
            return a;
        }
        S x1, y1;
        U d = extendedEuclidean(b, a % b, x1, y1);
        x = y1;
        y = x1 - y1 * (a / b);
        return d;
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
    template<typename T, std::enable_if_t<std::is_integral<T>::value || std::is_floating_point<T>::value, bool> = true>
    int sgn(T x) {
        return (x > T(0)) - (x < T(0));
    }

    // O(1)
    // Sets `a` to the minimum of `a` and `b`
    // @note `a` is provided by reference
    template<typename T, std::enable_if_t<std::is_integral<T>::value || std::is_floating_point<T>::value, bool> = true>void pMin(T &a, T b) {
        if (b < a) a = b;
    }

    // O(1)
    // Sets `a` to the maximum of `a` and `b`
    // @note `a` is provided by reference
    template<typename T, std::enable_if_t<std::is_integral<T>::value || std::is_floating_point<T>::value, bool> = true>
    void pMax(T &a, T b) {
        if (b > a) a = b;
    }

    /************************************************
     *                OTHER UTILITIES               *
     ************************************************/

    // @returns The number of characters of the printed representation of `obj`
    template<typename T> std::string repr(const T &obj) {
        std::stringstream s;
        s << obj;
        return s.str();
    }

    template<typename T> std::basic_string<T> repeat(const std::basic_string<T>& input, size_t num) {
        std::basic_stringstream<T> os;
        std::fill_n(std::ostream_iterator<std::basic_string<T>, T>(os), num, input);
        return os.str();
    }

    template<typename T> std::basic_string<T> operator*(std::basic_string<T> str, std::size_t n) {
        return repeat(std::move(str), n);
    }


    // std::wstring repeat(const std::wstring& input, size_t num) {
    //     std::wstringstream os;
    //     std::fill_n(std::ostream_iterator<std::wstring, wchar_t>(os), num, input);
    //     return os.str();
    // }

    // std::wstring operator*(std::wstring str, std::size_t n) {
    //     return repeat(std::move(str), n);
    // }

    // whats the units?
    double timeNow() {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
        ).count() / 1000.;
    }

    /************************************************
     *               ARRAY USEFUL IMPL              *
     ************************************************/

    // O(n)
    template<typename It>
    typename std::iterator_traits<It>::value_type sum(const It &begin, const It &end) {
        return accumulate(
            begin,
            end,
            typename std::iterator_traits<It>::value_type(0)
        );
    }

    // O(n)
    template<typename It1, typename It2>
    void make_prefix(const It1 &begin, const It1 &end, It2 out) {
        typename std::iterator_traits<It1>::value_type total = typename std::iterator_traits<It1>::value_type(0);
        for (It1 it = begin; it != end; it++, out++) {
            total += *it;
            *out = total;
        }
    }

    // O(n)
    // *(out+i) = *(begin+i) - *(begin+i-1) when i, i-1 is valid
    // NOTE: *out is left untouched
    template<typename It1, typename It2>
    void make_suffix(const It1 &begin, const It1 &end, It2 out) {
        for (It1 it = begin; it != end; ) {
            out++;
            It1 old = it;
            it++;
            *out = *it - *old;
        }
    }

    /************************************************
     *               VECTOR LINALG OPS              *
     ************************************************/

    // O(n)
    template<typename T>
    std::vector<T> operator+(std::vector<T> &a, std::vector<T> &b) {
        assert(a.size() == b.size());
        std::vector<T> out(a.size());
        for (size_t i = 0; i < a.size(); ++i) out[i] = a[i] + b[i];
        return out;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator+=(std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) a[i] += b[i];
        return a;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator+(std::vector<T> &b, const T &lambda) {
        std::vector<T> out(b.size());
        for (size_t i = 0; i < b.size(); ++i) out[i] = lambda + b[i];
        return out;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator+=(std::vector<T> &b, const T &lambda) {
        for (size_t i = 0; i < b.size(); ++i) b[i] += lambda;
        return b;
    }

    // O(n)
    template<typename T>
    std::vector<T> operator-(std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        std::vector<T> out(a.size());
        for (size_t i = 0; i < a.size(); ++i) out[i] = a[i] - b[i];
        return out;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator-=(std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) a[i] -= b[i];
        return a;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator-(std::vector<T> &v, const T &lambda) {
        std::vector<T> out(v.size());
        for (size_t i = 0; i < v.size(); ++i) out[i] = lambda - v[i];
        return out;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator-=(std::vector<T> &v, T &lambda) {
        for (size_t i = 0; i < v.size(); ++i) v[i] -= lambda;
        return v;
    }

    // O(n)
    template<typename T>
    std::vector<T> operator*(std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        std::vector<T> out(a.size());
        for (size_t i = 0; i < a.size(); ++i) out[i] = a[i] * b[i];
        return out;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator*=(std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) a[i] *= b[i];
        return a;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator*(std::vector<T> &v, const T &lambda) {
        std::vector<T> out(v.size());
        for (size_t i = 0; i < v.size(); ++i) out[i] = lambda * v[i];
        return out;
    }
    // O(n)
    template<typename T>
    std::vector<T> operator*=(std::vector<T> &v, const T &lambda) {
        for (size_t i = 0; i < v.size(); ++i) v[i] *= lambda;
        return v;
    }

    /************************************************
     *           VECTOR BOOLEAN COMPARISON          *
     ************************************************/

    // O(n)
    template<typename T>
    bool operator==(std::vector<T> &a, std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return false;
        }
        return true;
    }

    // O(n)
    template<typename T>
    bool operator!=(std::vector<T> &a, std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] != b[i]) return true;
        }
        return false;
    }

    /************************************************
     *                 VECTOR READ                  *
     ************************************************/

    template<typename T>
    std::istream& operator>>(std::istream &in, std::vector<T> &v) {
        std::copy_n(std::istream_iterator<T>(in), v.size(), v.begin());
        // for (T &val : v) in >> val;
        return in;
    }

    /************************************************
     *           SEQUENCE LEXIGRAPHIC ORDER         *
     ************************************************/

    // O(n)
    template<typename T>
    bool operator<=(const std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] < b[i]) return true;
            else if (a[i] > b[i]) return false;
        }
        return true;
    }

    // O(n)
    template<typename T>
    bool operator<(const std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] < b[i]) return true;
            else if (a[i] > b[i]) return false;
        }
        return false;
    }

    // O(n)
    template<typename T>
    bool operator>=(const std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > b[i]) return true;
            else if (a[i] < b[i]) return false;
        }
        return true;
    }

    // O(n)
    template<typename T>
    bool operator>(const std::vector<T> &a, const std::vector<T> &b) {
        assert(a.size() == b.size());
        for (size_t i = 0; i < a.size(); ++i) {
            if (a[i] > b[i]) return true;
            else if (a[i] < b[i]) return false;
        }
        return false;
    }
};

#endif
