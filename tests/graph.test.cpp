#include "../src/Graph.h"
#include <iostream>
#include <fstream>
using namespace std;
using namespace DS;

void testGraph() {
	ifstream cin{ "graphin.txt" };
	int N, M;
	cin >> N >> M;

	Graph<int> graph(N, M, cin, false, true);
	auto a = graph.greedyColouring();
	cout << ~graph;
}

int main() {
	testGraph();
}