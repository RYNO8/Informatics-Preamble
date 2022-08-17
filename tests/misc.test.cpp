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
using namespace std;
using namespace DS;

void testFib() {
	Matrix<ModInt> fib1 = Grid<ModInt>(vector<vector<ModInt>>{
		{ 1, 1 },
		{ 1, 0 }
	});
	for (int i = 0; i <= 50; ++i) cout << fib1.pow(i).getVal(0, 0) << "\n";
}

void testPopcount(ll trials = 1e9) {
	for (ll i = 0; i < trials; ++i) {
		uint32_t x = uint_dis(rng);
		assert(popcount(x) == __builtin_popcount(x));
		uint64_t y = ull_dis(rng);
		assert(popcount(y) == __builtin_popcountll(y));
	}
}

int main() {
	//testFib();
	testPopcount(1e4);
}