// alternatively `#include <bits/stdc++.h>`, but my compiler doesn't support that :'(
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include <iterator>
#include <utility>
#include <numeric>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <bitset>
#include <random>
#include <chrono>
#include <math.h>
#include <time.h>
#include <string>
#include <vector>
#include <limits>
#include <queue>
#include <stack>
#include <set>
#include <map>

using namespace std;
using namespace std::chrono;

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
constexpr ll MOD = 1000000007;

vector<pair<int, int>> DIRS_RECTILINEAR = { {0, 1}, {1, 0}, {0, -1}, {-1, 0} };
vector<pair<int, int>> DIRS_DIAG = { { 1, 1 }, { 1, -1 }, { -1, 1 }, { -1, -1 } };
vector<pair<int, int>> DIRS_ALL = { {0, 1}, {1, 0}, {0, -1}, {-1, 0}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1} };
vector<pair<int, int>> DIRS = DIRS_RECTILINEAR;

// Sets `a` to the minimum of `a` and `b`
// @note `a` is provided by reference
template<typename T> void pMin(T& a, T b) {
	if (b < a) a = b;
}

// Sets `a` to the maximum of `a` and `b`
// @note `a` is provided by reference
template<typename T> void pMax(T& a, T b) {
	if (b > a) a = b;
}

// Displays a pair
template<typename A, typename B> ostream& operator<<(ostream& os, const pair<A, B>& c) {
	os << "(" << c.first << ", " << c.second << ")";
	return os;
}

// Displays a vector
template<typename T> ostream& operator<<(ostream& os, const vector<T>& v) {
	os << "{ ";
	for (int i = 0; i < (int)v.size(); ++i) {
		os << v[i];
		if (i != v.size() - 1) os << ", ";
	}
	os << " }";
	return os;
}

// Displays an array
#define printArr0(a, N) for (int i = 0; i < N; ++i) cout << a[i] << " ";
#define printArr1(a, N) for (int i = 1; i <= N; ++i) cout << a[i] << " ";
#define printArr(a, N) printArr0(a, N);

// Displays a queue
template<typename T> ostream& operator<<(ostream& os, queue<T>& q) {
	os << "< ";
	for (int i = 0; i < (int)q.size(); ++i) {
		os << q.front();
		q.push(q.front());
		q.pop();
		if (i != q.size() - 1) os << ", ";
	}
	os << " >";
	return os;
}

// Displays a deque
template<typename T> ostream& operator<<(ostream& os, deque<T>& q) {
	os << "< ";
	for (int i = 0; i < (int)q.size(); ++i) {
		os << q.front();
		q.push_back(q.front());
		q.pop_front();
		if (i != q.size() - 1) os << ", ";
	}
	os << " >";
	return os;
}
