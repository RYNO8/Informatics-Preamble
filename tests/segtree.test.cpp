#include "../src/Segtree.h"
#include "../src/Ranges.h"
using namespace std;
using namespace DS;

// void testSegtree_manual() {
//     ll N;
//     cin >> N;

//     Segtree<ll> t(N);
//     cout << "Enter operation: set (S), add (A), min (m), max (M), total (T), quit (Q)\n\n";

//     while (true) {
//         char op;
//         cin >> op;

//         ll tl, tr, x;
//         switch (op) {
//         case 'S':
//             cin >> tl >> tr >> x;
//             t = *t.set(tl, tr, x);
//             break;
//         case 'A':
//             cin >> tl >> tr >> x;
//             t = *t.add(tl, tr, x);
//             break;
//         case 'm':
//             cin >> tl >> tr;
//             cout << t.getMin(tl, tr) << '\n';
//             break;
//         case 'M':
//             cin >> tl >> tr;
//             cout << t.getMax(tl, tr) << '\n';
//             break;
//         case 'T':
//             cin >> tl >> tr;
//             cout << t.query(tl, tr) << '\n';
//             break;
//         case 'Q':
//             return;
//         }

//         cout << t << '\n';

//         ll iMin = t.getMinIndexFirst();
//         assert(t[iMin] == t.getMin() && t.getMin(0, iMin - 1) > t[iMin]);
//         cout << string(int(iMin) * 2, ' ') << "^ min\n";

//         ll iMax = t.getMaxIndexFirst();
//         assert(t[iMax] == t.getMax() && t.getMax(0, iMax - 1) < t[iMax]);
//         cout << string(int(iMax) * 2, ' ') << "^ max\n";
//     }
// }

struct fakeSegtree {
    int N;
    vector<int> t;
    fakeSegtree(int _N, int initVal) : N(_N) {
        t = vector<int>(N, initVal);
    }

    void set(int l, int r, int x) {
        for (int i = l; i <= r; ++i) t[i] = x;
    }

    void add(int l, int r, int x) {
        for (int i = l; i <= r; ++i) t[i] += x;
    }

    int sum(int l, int r) {
        return accumulate(t.begin() + l, t.begin() + r + 1, 0);
    }

    int min(int l, int r) {
        return *min_element(t.begin() + l, t.begin() + r + 1);
    }

    int max(int l, int r) {
        return *max_element(t.begin() + l, t.begin() + r + 1);
    }

    int minIndex() {
        return min_element(t.begin(), t.end()) - t.begin();
    }

    int maxIndex() {
        return max_element(t.begin(), t.end()) - t.begin();
    }

    void extendRight(int x) {
        int s = t.size();
        for (int rep = 0; rep < s; ++rep) {
            t.push_back(x);
        }
    }

    size_t size() {
        return t.size();
    }
};

void testSegtree_auto() {
    for (int tI = 0; tI < 20; ++tI) {
        cout << "REP: " << tI << "\n";

        int N = 5;
        int initVal = 69;
        fakeSegtree t(N, initVal);
        ll shift = -rand();
        
        SumSegtree<ll> *tSum = new SumSegtree<ll>(0 + shift, N - 1 + shift, initVal);
        MaxSegtree<ll> *tMax = new MaxSegtree<ll>(0 + shift, N - 1 + shift, initVal);
        MinSegtree<ll> *tMin = new MinSegtree<ll>(0 + shift, N - 1 + shift, initVal);

        for (int rep = 0; rep < 10000; ++rep) {
            int opt = rand() % 8;
            int l = rand() % t.size(), r = rand() % t.size();
            if (l > r) swap(l, r);
            int x = (rand() % 100) - 50;

            if (opt == 0) {
                // set
                t.set(l, r, x);
                tSum->set({l + shift, r + shift}, x);
                tMax->set({l + shift, r + shift}, x);
                tMin->set({l + shift, r + shift}, x);

            } else if (opt == 1) {
                // add
                t.add(l, r, x);
                tSum->add({l + shift, r + shift}, x);
                tMax->add({l + shift, r + shift}, x);
                tMin->add({l + shift, r + shift}, x);

            } else if (opt == 2) {
                // get sum
                assert(t.sum(l, r) == tSum->query({l + shift, r + shift}));

            } else if (opt == 3) {
                // get min
                assert(t.min(l, r) == tMin->query({l + shift, r + shift}));

            } else if (opt == 4) {
                // get max
                assert(t.max(l, r) == tMax->query({l + shift, r + shift}));

            } else if (opt == 5) {
                // get index of min
                function<bool(MinSegtree<ll>*)> isMin = [&](MinSegtree<ll>* node) {
                    return node->query() == tMin->query();
                };
                assert(t.minIndex() == tMin->findFirst(isMin)->l() - shift);

            } else if (opt == 6) {
                // get index of max
                function<bool(MaxSegtree<ll>*)> isMax = [&](MaxSegtree<ll>* node) {
                    return node->query() == tMax->query();
                };
                assert(t.maxIndex() == tMax->findFirst(isMax)->l() - shift);

            } else if (opt == 7 && 2 * t.size() < 5e4) {
                // extend right
                t.extendRight(x);
                tSum = tSum->extendRight(x);
                tMax = tMax->extendRight(x);
                tMin = tMin->extendRight(x);
            }
        }
    }
}

signed main() {
    srand(time(NULL));
    //testSegtree_manual();
    testSegtree_auto();
}
