#include "../src/Constants.h"
#include "../src/ModInt.h"
using namespace std;

void testModInt() {
	auto a = ModInt(-20);
	cout << a << "\n";
	cout << a + ModInt(2) << "\n";
	cout << ModInt(2).pow(10) << "\n";
	cout << "\n";
}

signed main() {
	testModInt();
}