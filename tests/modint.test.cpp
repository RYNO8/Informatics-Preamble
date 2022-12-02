#include "../src/ModInt.h"
#include "../src/Constants.h"
#include "../src/Util.h"
#include <vector>
using namespace DS;
using namespace std;

void testModInt() {
    for (int rep = 0; rep < 1000000; ++rep) {
        auto a = ModInt<MOD>(ull_dis(rng));
        assert(a * a.modInv() == 1LL);
    }

    for (int rep = 0; rep < 1000000; ++rep) {
        auto n = ModInt<1000000007>(ull_dis(rng));
        auto d = ModInt<1000000007>(ull_dis(rng));
        vector<ModInt<1000000007>> ans = n/d;
        assert(ans.size() == gcd((uintmax_t)1000000007, (uintmax_t)d));
        for (auto val : ans) {
            assert(val * d == n);
        }
    }

    assert(ModInt<18>(9)/ModInt<18>(15) == vector<ModInt<18>>({15, 3, 9}));
    
    for (int rep = 0; rep < 1000000; ++rep) {
        auto a = ModInt<1000000007>(ull_dis(rng));
        assert(a.pow(1000000006) == 1LL);
    }

    assert(ModInt<6>(5) % 3 == 2);
}

void testModInt_bounds() {
    auto a = ModInt<numeric_limits<uint>::max()>(-1);
    assert(a == (ll)numeric_limits<uint>::max() - 1);
    assert((a * a) == 1LL);

    auto b = ModInt<numeric_limits<uint>::max()>(numeric_limits<uint>::max() - 1);
    assert(b == (ll)numeric_limits<uint>::max() - 1);
    assert((b * b) == 1LL);

    assert(ModInt<3>(numeric_limits<int>::min()) == 1LL);
    assert(ModInt<3>(numeric_limits<int>::max()) == 1LL);
    assert(ModInt<3>(numeric_limits<uint>::max()) == 0LL);
    assert(ModInt<3>(numeric_limits<ll>::min()) == 1LL);
    assert(ModInt<3>(numeric_limits<ll>::max()) == 1LL);
    assert(ModInt<3>(numeric_limits<ull>::max()) == 0LL);

    assert(ModInt<4294967295U>(numeric_limits<int>::min()) == 2147483647LL);
    assert(ModInt<4294967295U>(numeric_limits<int>::max()) == 2147483647LL);
    assert(ModInt<4294967295U>(numeric_limits<uint>::max()) == 0LL);
    assert(ModInt<4294967295U>(numeric_limits<ll>::min()) == 2147483647LL);
    assert(ModInt<4294967295U>(numeric_limits<ll>::max()) == 2147483647LL);
    assert(ModInt<4294967295U>(numeric_limits<ull>::max()) == 0LL);

    assert(ModInt<4294967294U>(numeric_limits<uint>::max()) == 1LL);
    assert(ModInt<4294967294U>(numeric_limits<ull>::max()) == 3LL);
}

signed main() {
    intmax_t largest = MOD - 1;
    assert((largest * largest) / largest == largest);

    assert(numeric_limits<uintmax_t>::max() > numeric_limits<intmax_t>::max());

    testModInt();
    testModInt_bounds();
}