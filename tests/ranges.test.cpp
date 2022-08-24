#include "../src/Constants.h"
#include "../src/Ranges.h"
using namespace std;
using namespace DS;

void testRanges() {
	/*vector<int> a = { 0,6,4, 4, 0, -1 };
	coordCompressTransform(a.begin(), a.end());
	cout << a << "\n";

	int b[] = { 0,6,4, 4, 0, -1 };
	coordCompressTransform(begin(b), end(b));
	printArr0(b, 6);
	cout << "\n";*/

	/*Range<int> r = {1, 2};
	cout << r << " " << (r < r);*/

	Ranges<int> r = { {1, 2}, {3, 4} };
	cout << r.size();
}


signed main() {
	testRanges();
}