#include "../src/Polynomial.h"
#include "../src/Constants.h"
using namespace DS;
using namespace std;

void testPolynomial() {
    Polynomial<ll> f({ -1, 1 });
    cout << f << "\n";
}

signed main() {
    testPolynomial();
}