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


int main() {
    Matrix<ModInt> fib1 = Grid<ModInt>(vector<vector<ModInt>>{
		{ 1, 1 },
		{ 1, 0 }
	});
	for (int i = 0; i <= 50; ++i) cout << fib1.pow(i).getVal(0, 0) << "\n";
}