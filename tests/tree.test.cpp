#include "../src/Tree.h"
using namespace std;
using namespace DS;


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

signed main() {
	testTree();
}