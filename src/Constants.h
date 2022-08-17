/************************************************
 *                    IMPORTS                   *
 ************************************************/

#include <bits/stdc++.h>
// #include <unordered_map>
// #include <unordered_set>
// #include <functional>
// #include <algorithm>
// #include <iostream>
// #include <assert.h>
// #include <iterator>
// #include <utility>
// #include <numeric>
// #include <cstddef>
// #include <fstream>
// #include <iomanip>
// #include <sstream>
// #include <climits>
// #include <bitset>
// #include <random>
// #include <chrono>
// #include <math.h>
// #include <time.h>
// #include <string>
// #include <vector>
// #include <limits>
// #include <queue>
// #include <stack>
// #include <set>
// #include <map>
#ifdef _MSC_VER
#	include <intrin.h>
#	define __builtin_popcount __popcnt
#	define __builtin_popcountll __popcntll
#endif

using namespace std;
using namespace std::chrono;

/************************************************
 *                   CONSTANTS                  *
 ************************************************/
typedef long double ld;
constexpr ld PI =    3.1415926535897932384626433832795;
constexpr ld E =     2.7182818284590452353602874713527;
constexpr ld PHI =   1.6180339887498948482045868343656;
constexpr ld GAMMA = 0.5772156649015328606065120900824;

typedef long long ll;
constexpr ll SIXTYNINE = 69;
constexpr ll MAXR = 1000;
constexpr ll MAXC = 1000;
constexpr ll MAXN = 1 << 19;
constexpr ll LOG_MAXN = 19;
constexpr ll SQRT_MAXN = 724; // or 300
constexpr ll MOD = 1000000007;
constexpr ll MOD_POW = 1; // highest value such that MOD-1 is divisible by 2^MOD_POW

typedef unsigned long long ull;

vector<pair<int, int>> DIRS_RECTILINEAR = { {0, 1}, {1, 0}, {0, -1}, {-1, 0} };
vector<pair<int, int>> DIRS_DIAG = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
vector<pair<int, int>> DIRS_ALL = { {0, 1}, {1, 0}, {0, -1}, {-1, 0}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1} };
vector<pair<int, int>> DIRS = DIRS_RECTILINEAR;

random_device rd;
mt19937_64 rng(rd());
uniform_int_distribution<int> int_dis(0, numeric_limits<int>::max());
uniform_int_distribution<ll> ll_dis(0, numeric_limits<ll>::max());
uniform_int_distribution<ull> ull_dis(0, numeric_limits<ull>::max());
uniform_real_distribution<ld> prob_dist(0, 1);

/************************************************
 *                    DISPLAY                   *
 ************************************************/
// Displays a pair:
// (A, B)
template<typename A, typename B> ostream& operator<<(ostream& out, const pair<A, B>& c) {
	out << "(" << c.first << ", " << c.second << ")";
	return out;
}

// Displays a vector:
// [ T, T, ..., T ]
template<typename T> ostream& operator<<(ostream& out, const vector<T>& v) {
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
#define printArr0(a, N) for (int i = 0; i < N; ++i) cout << a[i] << " ";
#define printArr1(a, N) for (int i = 1; i <= N; ++i) cout << a[i] << " ";
#define printArr(a, N) printArr0(a, N);

// Displays a queue:
// < T, T, ..., T >
template<typename T> ostream& operator<<(ostream& out, const queue<T>& q) {
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
template<typename T> ostream& operator<<(ostream& out, const stack<T>& s) {
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
template<typename T> ostream& operator<<(ostream& out, const deque<T>& q) {
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
template<typename A, typename B> ostream& operator<<(ostream& out, const map<A, B> &m) {
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
template<typename A, typename B> ostream& operator<<(ostream& out, const unordered_map<A, B>& m) {
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
template<typename T> ostream& operator<<(ostream& out, const set<T>& m) {
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
template<typename T> ostream & operator<<(ostream & out, const unordered_set<T> & m) {
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
template<typename T> ostream & operator<<(ostream & out, const multiset<T> & m) {
	out << "{ ";
	for (auto it = m.begin(); it != m.end(); ) {
		out << *it;
		if (++it != m.end()) out << ", ";
	}
	out << " }";
	return out;
}

/************************************************
 *                   UTILITIES                  *
 ************************************************/

// O(log a + log b)
// @note `a` and `b` should be non-negative
template<typename T> T gcd(T a, T b) {
	if (b == T(0)) return a;
	else return gcd(b, a % b);
}

// O(log a + log b)
template<typename T> T gcdSlow(T x, T y) {
	while (y != 0) {
		T temp = y;
		y = x % temp;
		x = temp;
	}
	return x;
}

// O(log a + log b)
// @note May overflow integer
template<typename T> T lcm(T a, T b) {
	return (a / gcd(a, b)) * b;
}


// O(12)
// @returns the number of set bits in the binary represetnation of `A`
int32_t popcount(int32_t x) {
    x -= ((x >> 1) & 0b01010101010101010101010101010101);
    x = ((x >> 2) & 0b00110011001100110011001100110011) + (x & 0b00110011001100110011001100110011);
    return ((x + (x >> 4) & 0b00001111000011110000111100001111) * 0x1010101) >> 24;
}

uint32_t popcount(uint32_t x) {
    x -= ((x >> 1) & 0b01010101010101010101010101010101);
    x = ((x >> 2) & 0b00110011001100110011001100110011) + (x & 0b00110011001100110011001100110011);
    return ((x + (x >> 4) & 0b00001111000011110000111100001111) * 0x1010101) >> 24;
}

int64_t popcount(int64_t x) {
    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    return (x * 0x0101010101010101ULL) >> 56;
}

uint64_t popcount(uint64_t x) {
    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    return (x * 0x0101010101010101ULL) >> 56;
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
	stringstream s;
	s << obj;
	return s.str().size();
}