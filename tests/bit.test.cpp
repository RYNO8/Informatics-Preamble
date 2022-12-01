#include "../src/BIT.h"
#include <assert.h>
#include <vector>
#include <random>
#include <time.h>
using namespace std;
using namespace DS;

void testBIT_metadata() {
    BIT<int, RANGE_QUERY_POINT_UPDATE, 10, 11> b;
    assert(b.size() == 110);
    assert(b.dimensions() == 2);
    assert(b.shape() == vector<size_t>({10, 11}));
}

void testBIT_rqpu() {
    // if it works in 2 dimensions, it probably works for 
    // all other dimensions
    BIT<int, RANGE_QUERY_POINT_UPDATE, 10, 11> b;
    int fakeBit[10][11] = {};

    for (int rep = 0; rep < 100000; ++rep) {
        int opt = rand() % 3;
        if (opt == 0) {
            // update index
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            int v = rand() % 69;
            fakeBit[i1][i2] += v;
            b.addIndex(i1, i2, v);
        } else if (opt == 1) {
            // query range
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
        } else {
            // query point
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            assert(fakeBit[i1][i2] == b.queryIndex(i1, i2));
        }
    }
}

void testBIT_pqru() {
    // if it works in 2 dimensions, it probably works for 
    // all other dimensions
    BIT<int, POINT_QUERY_RANGE_UPDATE, 10, 11> b;
    int fakeBit[10][11] = {};

    for (int rep = 0; rep < 100000; ++rep) {
        int opt = rand() % 3;
        if (opt == 0) {
            // query
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            
            assert(fakeBit[i1][i2] == b.queryIndex(i1, i2));
        } else if (opt == 1) {
            // update range
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
            b.addRange(l1, r1, l2, r2, v);
        } else if (opt == 2) {
            // update point
            int i1 = rand() % 10;
            int i2 = rand() % 11;
            int v = rand() % 69;
            fakeBit[i1][i2] += v;
            b.addIndex(i1, i2, v);
        }
    }
}

int main() {
    srand(time(NULL));
    testBIT_metadata();
    testBIT_rqpu();
    testBIT_pqru();
}