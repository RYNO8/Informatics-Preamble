#include "../src/Constants.h"
#include "../src/Polynomial.h"
using namespace std;

void testPolynomial() {
    Polynomial<ll> f({ -1, 1 });
    cout << f << "\n";
}

signed main() {
    testPolynomial();
}