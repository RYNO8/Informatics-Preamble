#include "../src/BIT.h"
#include "../src/CHT.h"
#include "../src/Constants.h"
#include "../src/Geometry.h"
#include "../src/Grid.h"
#include "../src/Graph.h"
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
    constexpr uint mod = 1e9;
    Matrix<ModInt<mod>> fib1 = Grid<ModInt<mod>>(vector<vector<ModInt<mod>>>{
        { 1, 1 },
        { 1, 0 }
    });
    assert(fib1.pow(100000).getVal(1, 0) == 428746875);
}


int main() {
    testFib();
}