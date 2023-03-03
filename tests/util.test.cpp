#include "../src/Util.h"
#include "../src/Constants.h"
using namespace DS;
using namespace std;

void testPopcount() {
    for (ll i = 0; i < 1e4; ++i) {
        uint32_t x = uint_dis(rng);
        assert(popcount(x) == __builtin_popcount(x));
        uint64_t y = ull_dis(rng);
        assert(popcount(y) == __builtin_popcountll(y));
    }
}

void testGCD_speed() {
    uint64_t ans = 0; // prevent compiler optimisation

    auto start1 = timeNow();
    for (ll i = 0; i < 1e8; ++i) {
        uint64_t x = ull_dis(rng);
        uint64_t y = ull_dis(rng);
        ans ^= gcd(x, y);
    }
    auto stop1 = timeNow();
    cout << stop1 - start1 << "\n";

    // auto start2 = timeNow();
    // for (ll i = 0; i < 1e8; ++i) {
    //     uint64_t x = ull_dis(rng);
    //     uint64_t y = ull_dis(rng);
    //     ans ^= gcdSlow(x, y);
    // }
    // auto stop2 = timeNow();
    // cout << stop2 - start2 << "\n";

    // auto start3 = timeNow();
    // for (ll i = 0; i < 1e8; ++i) {
    //     uint64_t x = ull_dis(rng);
    //     uint64_t y = ull_dis(rng);
    //     ans ^= __gcd(x, y);
    // }
    // auto stop3 = timeNow();
    // cout << stop3 - start3 << "\n";

    cout << ans << "\n";
}

void testVec() {
    vector<int> a = {1, 2, 3};
    vector<int> b = {4, 6, 1};
    assert(a <= b);
    a += b;
    assert(a == vector<int>({5, 8, 4}));
    assert(a - b == vector<int>({1, 2, 3}));
    assert(a != b);
}

int main() {
    testPopcount();
    //testGCD_speed();
    testVec();
}