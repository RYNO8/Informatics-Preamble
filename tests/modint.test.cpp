#include "../src/ModInt.h"
#include "../src/Constants.h"
#include "../src/Util.h"
#include <vector>
using namespace DS;
using namespace std;

void testModInt() {
    for (int rep = 0; rep < 1000000; ++rep) {
        auto a = ModInt<MOD>(ull_dis(rng));
        assert(a * a.inv() == 1LL);
        assert(a.pow(1000000006) == 1LL);
        assert(a.inv() == a.pow(-1));
    }

    for (int rep = 0; rep < 1000000; ++rep) {
        auto n = ModInt<1000000007>(ull_dis(rng));
        auto d = ModInt<1000000007>(ull_dis(rng));
        vector<ModInt<1000000007>> ans = n.allCongreunt(d);
        assert(ans.size() == gcd((uintmax_t)1000000007, d.num()));
        for (auto val : ans) {
            assert(val * d == n);
        }
    }
    assert(ModInt<18>(9).allCongreunt(ModInt<18>(15)) == vector<ModInt<18>>({15, 3, 9}));

    assert(repr(ModInt<4206969>(69)) == "69 ( mod 4206969 )");
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

void testModInt_atcoder() {
    ModInt<11> a = 10;
    ModInt<11> b(3);

    // equal
    assert(a == 21);
    assert(a == -1);
    assert(-1 == a);

    // negative
    assert(-b == 8);
    assert(0 - b == 8);

    // plus
    assert(a + b == 2);  // (10 + 3) mod 11
    assert(1 + a == 0);

    // minus
    assert(a - b == 7);  // (10 - 3) mod 11
    assert(b - a == 4);

    // mul
    assert(a * b == 8);  // (10 * 3) mod 11

    // inv
    assert(b.inv() == 4);  // (3 * 4) mod 11 == 1

    // div
    assert(a / b == 7);  // (10 * 4) mod 11

    // +=, -=, *=
    a += b;
    assert(a == 2 && b == 3);
    a -= b;
    assert(a == 10 && b == 3);
    a *= b;
    assert(a == 8 && b == 3);

    // pow
    assert(ModInt<11>(2).pow(4) == 5);  // 16 mod 11

    // get mod
    assert(ModInt<11>::mod() == 11 && a.mod() == 11);
}

signed main() {
    assert(is_my_integral<ModInt<10>>::value);

    intmax_t largest = MOD - 1;
    assert((largest * largest) / largest == largest);

    assert(numeric_limits<uintmax_t>::max() > numeric_limits<intmax_t>::max());

    testModInt();
    testModInt_bounds();
    testModInt_atcoder();
}