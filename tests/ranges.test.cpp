#include "../src/Constants.h"
#include "../src/Ranges.h"
using namespace std;

void testRanges() {
	vector<int> a = { 0,6,4, 4, 0, -1 };
	coordCompressTransform(a);
	cout << a;
}


signed main() {
	testRanges();
}