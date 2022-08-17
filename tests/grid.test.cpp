#include "../src/Grid.h"
#include <iostream>
#include <fstream>
using namespace std;
using namespace DS;

void testGrid() {
	ifstream cin{ "gridin.txt" };

	int R, C;
	cin >> R >> C;


	Grid<ll> g(R, C, cin);
	cout << g;
	auto b = g.shortestPath(0, 0, 0, 2, g.compLE(4), DIRS_RECTILINEAR);
}

signed main() {
	testGrid();
}