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


void testSegtree_auto() {
    for (int t = 0; t < 20; ++t) {
        cout << "REP: " << t << "\n";

        int fakeSegtree[10000] = {};
        ll shift = -rand();
        fill(begin(fakeSegtree), end(fakeSegtree), 69);
        SumSegtree<ll> tSum(0 + shift, 9999 + shift, 69);
        MaxSegtree<ll> tMax(0 + shift, 9999 + shift, 69);
        MinSegtree<ll> tMin(0 + shift, 9999 + shift, 69);

        for (int rep = 0; rep < 10000; ++rep) {
            int opt = rand() % 7;
            int l = rand() % 10000, r = rand() % 10000;
            if (l > r) swap(l, r);

            if (opt == 0) {
                // set
                int x = (rand() % 100) - 50;
                for (int i = l; i <= r; ++i) fakeSegtree[i] = x;
                tSum.set({l + shift, r + shift}, x);
                tMax.set({l + shift, r + shift}, x);
                tMin.set({l + shift, r + shift}, x);
            } else if (opt == 1) {
                // add
                int x = (rand() % 100) - 50;
                //for (int i = l; i <= r; ++i) fakeSegtree[i] += x;
                tSum.add({l + shift, r + shift}, x);
                tMax.add({l + shift, r + shift}, x);
                tMin.add({l + shift, r + shift}, x);
            } else if (opt == 2) {
                // get min
                assert(*min_element(fakeSegtree + l, fakeSegtree + r + 1) == tMin.query({l + shift, r + shift}));
            } else if (opt == 3) {
                // get max
                assert(*max_element(fakeSegtree + l, fakeSegtree + r + 1) == tMax.query({l + shift, r + shift}));
            } else if (opt == 4) {
                // get index of min
                function<bool(MinSegtree<ll>*)> isMin = [&](MinSegtree<ll>* node) {
                    return node->query() == tMin.query();
                };
                assert(min_element(begin(fakeSegtree), end(fakeSegtree)) - begin(fakeSegtree) == tMin.findFirst(isMin)->l() - shift);
            } else if (opt == 5) {
                // get index of max
                function<bool(MaxSegtree<ll>*)> isMax = [&](MaxSegtree<ll>* node) {
                    return node->query() == tMax.query();
                };
                assert(max_element(begin(fakeSegtree), end(fakeSegtree)) - begin(fakeSegtree) == tMax.findFirst(isMax)->l() - shift);
            } else if (opt == 6) {
                int expected = 0;
                for (int i = l; i <= r; ++i) expected += fakeSegtree[i];
                assert(expected == tSum.query({l + shift, r + shift}));
            }
        }
    }
}

signed main() {
    srand(time(NULL));
    //testSegtree_manual();
    testSegtree_auto();
}
