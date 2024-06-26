#ifndef CONSTANTS_H
#define CONSTANTS_H

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
#include <intrin.h>
#define __builtin_popcount __popcnt
#define __builtin_popcountll __popcntll
#endif
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>

/************************************************
 *                   CONSTANTS                  *
 ************************************************/

namespace DS {
template<class T, class Compare = std::less<T>>
using ordered_set = __gnu_pbds::tree<
    T, __gnu_pbds::null_type, Compare, __gnu_pbds::rb_tree_tag,
    __gnu_pbds::tree_order_statistics_node_update>;

typedef uint32_t uint;

typedef long double ld;
constexpr ld PI    = 3.1415926535897932384626433832795;
constexpr ld E     = 2.7182818284590452353602874713527;
constexpr ld PHI   = 1.6180339887498948482045868343656;
constexpr ld GAMMA = 0.5772156649015328606065120900824;

typedef int64_t ll;
constexpr ll SIXTYNINE = 69;
constexpr ll MAXR      = 1000;
constexpr ll MAXC      = 1000;
constexpr ll MAXN      = 1 << 19;
constexpr ll LOG_MAXN  = 19;
constexpr ll SQRT_MAXN = 724; // or 300
constexpr ll MOD       = 1000000007;
constexpr ll MOD_POW   = 1; // highest value such that MOD-1 is divisible by 2^MOD_POW

typedef uint64_t ull;

std::pair<int, int> DIRS_RECTILINEAR[4] = {{0, 1}, {1, 0}, {0, -1}, {-1, 0}};
std::pair<int, int> DIRS_DIAG[4]        = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
std::pair<int, int> DIRS_ALL[8]         = {{0, 1}, {1, 0},  {0, -1}, {-1, 0},
                                           {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
std::pair<int, int> DIRS_KNIGHT[8]      = {{1, 2},  {2, 1},  {1, -2},  {-2, 1},
                                           {-1, 2}, {2, -1}, {-1, -2}, {-2, -1}};

std::random_device rd;
std::mt19937_64 rng(rd());
std::uniform_int_distribution<int> int_dis(0, std::numeric_limits<int>::max());
std::uniform_int_distribution<uint> uint_dis(0, std::numeric_limits<uint>::max());
std::uniform_int_distribution<ll> ll_dis(0, std::numeric_limits<ll>::max());
std::uniform_int_distribution<ull> ull_dis(0, std::numeric_limits<ull>::max());
std::uniform_real_distribution<ld> prob_dist(0, 1);

class MyIntegralIndicator {};

template<typename T>
using is_my_integral = typename std::conditional_t<
    std::is_integral<T>::value || std::is_base_of<MyIntegralIndicator, T>::value, std::true_type,
    std::false_type>;
template<typename T>
constexpr bool is_my_integral_v = is_my_integral<T>::value;

template<typename T>
using is_signed_int = typename std::conditional_t<
    std::is_integral<T>::value && std::is_signed<T>::value, std::true_type, std::false_type>;
template<typename T>
constexpr bool is_signed_int_v = is_signed_int<T>::value;

template<typename T>
using is_unsigned_int = typename std::conditional_t<
    std::is_integral<T>::value && std::is_unsigned<T>::value, std::true_type, std::false_type>;
template<typename T>
constexpr bool is_unsigned_int_v = is_unsigned_int<T>::value;

}; // namespace DS

#endif
