#include "../src/BIT.h"
#include "../src/CHT.h"
#include "../src/Constants.h"
#include "../src/Geometry.h"
#include "../src/Graph.h"
#include "../src/Grid.h"
#include "../src/Matrix.h"
#include "../src/ModInt.h"
#include "../src/Polynomial.h"
#include "../src/Ranges.h"
#include "../src/Segtree.h"
#include "../src/SqrtDecomp.h"
#include "../src/Tree.h"
#include "../src/Util.h"
using namespace std;
using namespace DS;

void testFib() {
    using M        = ModInt<(ll)1e9>;
    Matrix<M> fib1 = Grid<M>(vector<vector<M>>{{1, 1}, {1, 0}});
    assert(fib1.pow(100000).getVal(1, 0) == 428746875);
}

void testPrint() {
    ll a[] = {1, 2, 3, 4, 5};
    assert(repr(capture<ll *>{begin(a), end(a)}) == "[ 1 2 3 4 5 ]");
    vector<ll> v = vector<ll>({1, 2, 3, 4, 5});
    assert(repr(v) == "[ 1 2 3 4 5 ]");
}

int main() {
    testFib();
    testPrint();
}
