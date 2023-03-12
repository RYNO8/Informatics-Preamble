#include "../src/BIT.h"
#include "../src/Util.h"
#include <assert.h>
#include <vector>
#include <random>
#include <time.h>
using namespace std;
using namespace DS;

void testBIT_metadata() {
    BIT_RQPU<int, 2, 3> b;
    b.addIndex(0, 0, 0);
    b.addIndex(1, 0, 1);
    b.addIndex(2, 0, 2);
    b.addIndex(3, 1, 0);
    b.addIndex(4, 1, 1);
    b.addIndex(5, 1, 2);
    assert(b.size() == 6);
    assert(b.dimensions() == 2);
    assert(b.shape() == vector<size_t>({2, 3}));
    assert(b.index() == vector<vector<size_t>>({
        {0, 0},
        {0, 1},
        {0, 2},
        {1, 0},
        {1, 1},
        {1, 2},
    }));
    assert(b.data() == vector<int>({0, 1, 2, 3, 4, 5}));

    BIT_PQRU<int, 2> c;
    c.addIndex(0, 0);
    c.addIndex(1, 1);
    assert(c.data() == vector<int>({0, 1}));
}

void testBIT_rqpu() {
    // if it works in 2 dimensions, it probably works for 
    // all other dimensions
    BIT_RQPU<int, 10, 11> b;
    int fakeBit[10][11] = {};

    for (int rep = 0; rep < 100000; ++rep) {
        int opt = rand() % 4;
        if (opt == 0) {
            // range query
            int l1 = rand() % 10, r1 = rand() % 10;
            if (l1 > r1) swap(l1, r1);
            int l2 = rand() % 11, r2 = rand() % 11;
            if (l2 > r2) swap(l2, r2);
            int expectedTotal = 0;
            for (int i1 = l1; i1 <= r1; ++i1) {
                for (int i2 = l2; i2 <= r2; ++i2) {
                    expectedTotal += fakeBit[i1][i2];
                }
            }
            assert(expectedTotal == b.querySum(l1, r1, l2, r2));
        } else if (opt == 1) {
            // point query
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            assert(fakeBit[i1][i2] == b.queryIndex(i1, i2));
        } else if (opt == 2) {
            // point add-update
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            int v = rand() % 69;
            fakeBit[i1][i2] += v;
            b.addIndex(v, i1, i2);
        } else if (opt == 3) {
            // point set-update
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            int v = rand() % 69;
            fakeBit[i1][i2] = v;
            b.setIndex(v, i1, i2);
        }
    }
}

void testBIT_pqru() {
    // if it works in 2 dimensions, it probably works for 
    // all other dimensions
    BIT_PQRU<int, 10, 11> b;
    int fakeBit[10][11] = {};

    for (int rep = 0; rep < 100000; ++rep) {
        int opt = rand() % 4;
        if (opt == 0) {
            // point query
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            
            assert(fakeBit[i1][i2] == b.queryIndex(i1, i2));
        } else if (opt == 1) {
            // range add-update
            int l1 = rand() % 10, r1 = rand() % 10;
            if (l1 > r1) swap(l1, r1);
            int l2 = rand() % 11, r2 = rand() % 11;
            if (l2 > r2) swap(l2, r2);
            int v = rand() % 69;
            for (int i1 = l1; i1 <= r1; ++i1) {
                for (int i2 = l2; i2 <= r2; ++i2) {
                    fakeBit[i1][i2] += v;
                }
            }
            b.addRange(v, l1, r1, l2, r2);
        } else if (opt == 2) {
            // point add-update
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            int v = rand() % 69;
            fakeBit[i1][i2] += v;
            b.addIndex(v, i1, i2);
        } else if (opt == 3) {
            // point set-update
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            int v = rand() % 69;
            fakeBit[i1][i2] = v;
            b.setIndex(v, i1, i2);
        }
    }
}


void testBIT_display() {
    BIT_PQRU<int, 2, 2> b;
    b.addIndex(69, 0, 1);
    assert(repr(b) == "[ [ 0 69 ] [ 0 0 ] ]");
}

int main() {
    srand(time(NULL));
    testBIT_metadata();
    testBIT_rqpu();
    testBIT_pqru();
    testBIT_display();
}