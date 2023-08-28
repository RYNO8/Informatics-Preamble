#include "../src/Util.h"
#include "../src/Constants.h"
using namespace DS;
using namespace std;

void testDisplay() {
    assert(repr(vector<int>({})) == "[ ]");
    assert(repr(vector<int>({ 1, 2 })) == "[ 1 2 ]");

    assert(repr(array<int, 0>()) == "[ ]");
    assert(repr(array<int, 2>({ 1, 2 })) == "[ 1 2 ]");

    queue<int> q;
    assert(repr(q) == "< >");
    q.push(1);
    q.push(2);
    assert(repr(q) == "< 1, 2 >");

    stack<int> s;
    assert(repr(s) == "< >");
    s.push(2);
    s.push(1);
    assert(repr(s) == "< 1, 2 >");

    deque<int> dq;
    assert(repr(dq) == "< >");
    dq.push_back(1);
    dq.push_back(2);
    assert(repr(dq) == "< 1, 2 >");

    assert(repr(map<int, int>({})) == "{ }");
    assert(repr(map<int, int>({ { 1, 2 } })) == "{ 1: 2 }");
    assert(repr(unordered_map<int, int>({})) == "{ }");
    assert(repr(unordered_map<int, int>({ { 1, 2 } })) == "{ 1: 2 }");
    assert(repr(multimap<int, int>({})) == "{ }");
    assert(repr(multimap<int, int>({ { 1, 2 } })) == "{ 1: 2 }");
    assert(repr(unordered_multimap<int, int>({})) == "{ }");
    assert(repr(unordered_multimap<int, int>({ { 1, 2 } })) == "{ 1: 2 }");

    assert(repr(set<int>()) == "{ }");
    assert(repr(set<int>({ 1, 2 })) == "{ 1, 2 }");
    assert(repr(unordered_set<int>()) == "{ }");
    assert(
        repr(unordered_set<int>({ 1, 2 })) == "{ 1, 2 }" ||
        repr(unordered_set<int>({ 1, 2 })) == "{ 2, 1 }"
    );
    assert(repr(multiset<int>()) == "{ }");
    assert(repr(multiset<int>({ 1, 2 })) == "{ 1, 2 }");
    assert(repr(unordered_multiset<int>()) == "{ }");
    assert(
        repr(unordered_multiset<int>({ 1, 2 })) == "{ 1, 2 }" ||
        repr(unordered_multiset<int>({ 1, 2 })) == "{ 2, 1 }"
    );
}

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

    auto start1  = timeNow();
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
    vector<int> a = { 1, 2, 3 };
    vector<int> b = { 4, 6, 1 };
    assert(a <= b);
    a += b;
    assert(a == vector<int>({ 5, 8, 4 }));
    assert(a - b == vector<int>({ 1, 2, 3 }));
    assert(a != b);

    stringstream ss;
    ss << "3 5 0";
    vector<int> c(3);
    ss >> c;
    assert(c + 1 == b);
}

int main() {
    testDisplay();
    testPopcount();
    // testGCD_speed();
    testVec();
}
