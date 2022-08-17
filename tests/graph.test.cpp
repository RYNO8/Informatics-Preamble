// constants come first
#include "../src/Constants.h"
#include "../src/Graph.h"
using namespace std;

void testGraph() {
	ifstream cin{ "graphin.txt" };
	int N, M;
	cin >> N >> M;

	Graph<int> graph(N, M, cin, false, true);
	auto a = graph.greedyColouring();
	cout << ~graph;
}


signed main() {
	testGraph();
}