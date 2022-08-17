#include "../src/Constants.h"
#include "../src/Grid.h"
using namespace std;

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