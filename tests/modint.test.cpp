#include "../src/ModInt.h"
#include "../src/Constants.h"
using namespace DS;
using namespace std;

void testModInt() {
    for (int rep = 0; rep < 1000000; ++rep) {
        auto a = ModInt<>(ull_dis(rng));
        assert(a * a.modInv() == 1ll);
    }
}

void testModInt_bounds() {
    auto a = ModInt<numeric_limits<uint>::max()>(-1);
    assert(a == (ll)numeric_limits<uint>::max() - 1);
    assert((a * a) == 1ll);

    auto b = ModInt<numeric_limits<uint>::max()>(numeric_limits<uint>::max() - 1);
    assert(b == (ll)numeric_limits<uint>::max() - 1);
    assert((b * b) == 1ll);

    assert(ModInt<3>(numeric_limits<int>::min()) == 1ll);
    assert(ModInt<3>(numeric_limits<int>::max()) == 1ll);
    assert(ModInt<3>(numeric_limits<uint>::max()) == 0ll);
    assert(ModInt<3>(numeric_limits<ll>::min()) == 1ll);
    assert(ModInt<3>(numeric_limits<ll>::max()) == 1ll);
    assert(ModInt<3>(numeric_limits<ull>::max()) == 0ll);
}

signed main() {
    testModInt();
    testModInt_bounds();
}