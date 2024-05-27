#include "../src/Ranges.h"

#include "../src/Constants.h"
#include "../src/Util.h"

using namespace std;
using namespace DS;

void testCompression() {
    vector<int> a = { 0, 6, 4, 4, 0, -1 };
    coordCompressTransform(a.begin(), a.end());
    assert(a == vector<int>({ 1, 3, 2, 2, 1, 0 }));

    int b[] = { 0, 6, 4, 4, 0, -1 };
    coordCompressTransform(begin(b), end(b));
    int expectedB[] = { 1, 3, 2, 2, 1, 0 };
    assert(equal(begin(b), end(b), begin(expectedB)));

    array<int, 6> c = { 0, 6, 4, 4, 0, -1 };
    map<int, int> expectedC = { { -1, 0 }, { 0, 1 }, { 4, 2 }, { 6, 3 } };
    assert(coordCompressMap(c.begin(), c.end()) == expectedC);

    set<int> d = { -1, 0, 4, 6 };
    map<int, int> expectedD = { { -1, 0 }, { 0, 1 }, { 4, 2 }, { 6, 3 } };
    assert(coordCompressMap(d.begin(), d.end()) == expectedD);
}

void testRange() {
    // some methods (i.e. length, midpoint, covers, isCovered) are tested in segtree
    Range<int> a = { 1, 10 };
    Range<int> b = { 4, 5 };
    assert(a.covers(b));
    assert(b.coveredBy(a));
}

void testRanges() {
    Ranges<int> e = { { INT_MIN, INT_MAX } };
    Ranges<int> r1 = { { 1, 2 }, { 6, 9 } };
    Ranges<int> r1_not = { { INT_MIN, 0 }, { 3, 5 }, { 10, INT_MAX } };
    Ranges<int> r2 = { { -1, 1 }, { 5, 7 }, { 9, 20 }, { 40, 69 } };
    Ranges<int> r1_and_r2 = { { 1, 1 }, { 6, 7 }, { 9, 9 } };
    Ranges<int> r1_or_r2 = { { -1, 2 }, { 5, 20 }, { 40, 69 } };

    assert(r1.size() == 2);
    assert(r1.length() == 6);
    assert(!r1.covers(0) && r1.covers(1) && r1.covers(2) && !r1.covers(3) && !r1.covers(69));
    assert(r1.covers({ 7, 8 }) && r1.covers({ 7, 9 }) && r1.covers({ 2, 2 }));
    assert(!r1.covers({ 7, 10 }) && !r1.covers({ 3, 4 }) && !r1.covers({ 2, 3 }));

    assert(
        (r1 & r2) == r1_and_r2 && (r1 & r1_and_r2) == r1_and_r2 && (r2 & r1_and_r2) == r1_and_r2 &&
        (r1 & r1_not).empty()
    );

    assert((r1 | r2) == r1_or_r2 && (r1 | r1_not) == e);

    assert(~r1 == r1_not);
}

signed main() {
    testCompression();
    testRange();
    testRanges();
}
