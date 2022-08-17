// constants come first
#include "Constants.h"

// include the data structures you want
// TODO: make sure all orderings work (https://www.learncpp.com/cpp-tutorial/header-guards/)
#include "CHT.h"
#include "Geometry.h"
#include "Graph.h"
#include "Grid.h"
#include "Matrix.h" // must come after grid
#include "ModInt.h"
#include "Polynomial.h"
#include "Ranges.h"
#include "Segtree.h"
#include "SqrtDecomp.h"
#include "Tree.h" // must come after graph
using namespace std;


int main() {
    Matrix<ModInt> fib1 = Grid<ModInt>(vector<vector<ModInt>>{
		{ 1, 1 },
		{ 1, 0 }
	});
	for (int i = 0; i <= 50; ++i) cout << fib1.pow(i).getVal(0, 0) << "\n";
}