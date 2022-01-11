/* RYAN'S PREAMBLE
* This implements most common functions for grids, graphs and segtrees
* There are also some helpful templates to help with standard formats of questions, such as Google Kickstart/Codejam
* It is probably not the fastest or most memory efficient, because I have prioritised verbosity over constant factor optimisations
* I hope that once its completed, its size would still be less than the maximum file size allowed
* TODO: implement coordinate geometry, eww
* TODO: wrapper for sqrt decomp
* TODO: n dimensional fenwick tree
* TODO: dynamic tree for Tree.h, how??
* TODO: fast factorisation, fast isPrime, how?? Pollard rho was dodgy
* TODO: check through everything for consistency
* TODO: change arguments to references if possible
*/

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

void testGeometry() {
	Line<ll> a(Point<ll>(1, 1), Point<ll>(-1, -1));
	Line<ll> b(Point<ll>(1, -1), Point<ll>(-1, -1));
	cout << a.intersects(b);
}

void testGraph() {
	ifstream cin{ "graphin.txt" };
	int N, M;
	cin >> N >> M;

	Graph<int> graph(N, M, cin, false, true);
	auto a = graph.greedyColouring();
	cout << ~graph;
}

void testGrid() {
	ifstream cin{ "gridin.txt" };

	int R, C;
	cin >> R >> C;


	Grid<ll> g(R, C, cin);
	cout << g;
	auto b = g.shortestPath(0, 0, 0, 2, g.compLE(4), DIRS_RECTILINEAR);
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

void testMatrix() {
	Matrix<ModInt> fib1 = vector<vector<ModInt>>{
		{ 1, 1 },
		{ 1, 0 }
	};
	for (int i = 0; i <= 50; ++i) cout << fib1.pow(i).getVal(0, 0) << "\n";
}

void testModInt() {
	auto a = ModInt(-20);
	cout << a << "\n";
	cout << a + ModInt(2) << "\n";
	cout << ModInt(2).pow(10) << "\n";
	cout << "\n";
}

void testPolynomial() {
	Polynomial<ll> f({ -1, 1 });
	cout << f;
}

void testRanges() {
	vector<int> a = { 0,6,4, 4, 0, -1 };
	coordCompressTransform(a);
	cout << a;
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

		ll iMin = t.getMinIndexFirst();
		assert(t[iMin] == t.getMin() && t.getMin(0, iMin - 1) > t[iMin]);
		cout << string(int(iMin) * 2, ' ') << "^ min\n";

		ll iMax = t.getMaxIndexFirst();
		assert(t[iMax] == t.getMax() && t.getMax(0, iMax - 1) < t[iMax]);
		cout << string(int(iMax) * 2, ' ') << "^ max\n";
	}
}

void testSqrtDecomp() {
	int N, Q, bad = 0, arr[100005], curr[100005];
	vector<pair<int, int>> queries;

	cin >> N >> Q;
	for (int i = 1; i <= N; ++i) cin >> arr[i];
	for (int a, b, i = 0; i < Q; ++i) {
		cin >> a >> b;
		queries.push_back({ a, b });
	}

	coordCompressTransform(arr + 1, arr + N + 1);

	function<void(int)> add = [&](int i) {
		if (curr[arr[i]] % 2 == 0) ++bad;
		else --bad;
		curr[arr[i]]++;
	};

	function<void(int)> rem = [&](int i) {
		if (curr[arr[i]] % 2 == 0) ++bad;
		else --bad;
		curr[arr[i]]--;
	};

	function<bool(int, int)> answer = [&](int l, int r) {
		return bad == ((r - l + 1) % 2 == 1);
	};

	for (bool ans : chunkQueries(queries, add, rem, answer)) cout << (ans ? "YES" : "NO") << '\n';
}

signed main() {
	cin.tie(0); ios::sync_with_stdio(0);

	//testGeometry();
	//testGraph();
	//testGrid();
	//testTree();
	//testMatrix();
	//testModInt();
	//testPolynomial();
	//testRanges();
	//testSegtree();
	//testSqrtDecomp();
	set<int> a = { 1,3, 2 };
	cout << &a;
}