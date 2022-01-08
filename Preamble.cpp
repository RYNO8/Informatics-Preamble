/* RYAN'S PREAMBLE
* This implements most common functions for grids, graphs and segtrees
* There are also some helpful templates to help with standard formats of questions, such as Google Kickstart/Codejam
* It is probably not the fastest, because I have prioritised verbosity over constant factor optimisations
* I hope that once its completed, its size would still be less than the maximum file size allowed
* TODO: use iterators in Graph, so that it can be used through std algorithm
* e.g. find(g.begin(), g.end(), 1);
*/

// constants come first
#include "Constants.h"

// include the data structures you want
// TODO: make sure all orderings work (https://www.learncpp.com/cpp-tutorial/header-guards/)
#include "Grid.h"
#include "Graph.h"
#include "Tree.h" // must come after graph
#include "NumTheory.h" // must come after grid
#include "Segtree.h"
#include "CHT.h"

// not necessary
using namespace std;

void testGrid() {
	ifstream cin{ "gridin.txt" };

	int R, C;
	cin >> R >> C;


	Grid<ll> g(R, C, cin);
	cout << g;
	auto b = g.shortestPath(0, 0, 0, 2, g.compLE(4), DIRS_RECTILINEAR);
}

void testGraph() {
	ifstream cin{ "graphin.txt" };
	int N, M;
	cin >> N >> M;

	Graph<int> graph(N, M, cin, false, true);
	auto a = graph.greedyColouring();

}


void testTree() {
	ifstream fin{ "treein.txt" };
	int N;
	fin >> N;
	Graph<int> g(N, N - 1, fin, true, true);
	Tree<int> t(g);
	t.addNode(11);
	t.pprint();
	auto v = t.getPath(12, 6);
}

void testNumTheory() {
	auto a = ModInt(-20, 10);
	cout << a << "\n";
	cout << a + ModInt(2, 10) << "\n";
	cout << ModInt(2, 20).pow(10) << "\n";
	cout << "\n";

	Matrix<ModInt> fib1 = vector<vector<ModInt>>{
		{ 1, 1 },
		{ 1, 0 }
	};
	for (int i = 0; i <= 50; ++i) cout << fib1.pow(i).getVal(0, 0) << "\n";
}

void testSegtree() {
	ll N;
	cin >> N;

	Segtree<ll> t(N);
	cout << "Enter operation: set (S), add (A), min (m), max (M), total (T), quit (Q)\n\n";

	while (true) {
		char op;
		cin >> op;

		ll tl, tr, x;
		switch (op) {
		case 'S':
			cin >> tl >> tr >> x;
			t = *t.set(tl, tr, x);
			break;
		case 'A':
			cin >> tl >> tr >> x;
			t = *t.add(tl, tr, x);
			break;
		case 'm':
			cin >> tl >> tr;
			cout << t.getMin(tl, tr) << '\n';
			break;
		case 'M':
			cin >> tl >> tr;
			cout << t.getMax(tl, tr) << '\n';
			break;
		case 'T':
			cin >> tl >> tr;
			cout << t.getSum(tl, tr) << '\n';
			break;
		case 'Q':
			return;
		}

		cout << t << '\n';

		ll iMin = t.getMinIndex();
		assert(t[iMin] == t.getMin() && t.getMin(0, iMin - 1) > t[iMin]);
		cout << string(iMin * 2, ' ') << "^ min\n";

		ll iMax = t.getMaxIndex();
		assert(t[iMax] == t.getMax() && t.getMax(0, iMax - 1) < t[iMax]);
		cout << string(iMax * 2, ' ') << "^ max\n";
	}
}



signed main() {
	cin.tie(0); ios::sync_with_stdio(0);

	//testGrid();
	//testGraph();
	//testTree();
	//testNumTheory();
	testSegtree();

	//int a[] = { 1, 2, 3 };
	//printArr0(a, 3);
}