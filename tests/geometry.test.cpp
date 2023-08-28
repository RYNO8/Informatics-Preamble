// constants come first
#include "../src/Geometry.h"
#include "../src/Constants.h"
using namespace std;

void testGeometry() {
    Line<ll> a(Point<ll>(1, 1), Point<ll>(-1, -1));
    Line<ll> b(Point<ll>(1, -1), Point<ll>(-1, -1));
}

signed main() { testGeometry(); }
