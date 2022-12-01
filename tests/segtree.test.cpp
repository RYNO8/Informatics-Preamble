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
    int fakeSegtree[10000] = {};
    ll shift = 0; //-rand();
    fill(begin(fakeSegtree), end(fakeSegtree), 69);
    SumSegtree<ll> tSum(0 + shift, 10000 + shift, 69);
    MaxSegtree<ll> tMax(0 + shift, 10000 + shift, 69);
    MinSegtree<ll> tMin(0 + shift, 10000 + shift, 69);

    for (int rep = 0; rep < 10000; ++rep) {
        int opt = rand() % 7;
        int l = rand() % 10000, r = rand() % 10000;
        if (l > r) swap(l, r);

        if (opt == 0) {
            // set
            int x = (rand() % 100) - 50;
            for (int i = l; i <= r; ++i) fakeSegtree[i] = x;
            tSum.set(Range<ll>({l, r}) + shift, x);
            tMax.set(Range<ll>({l, r}) + shift, x);
            tMin.set(Range<ll>({l, r}) + shift, x);
        } else if (opt == 1) {
            // add
            int x = (rand() % 100) - 50;
            for (int i = l; i <= r; ++i) fakeSegtree[i] += x;
            tSum.add(Range<ll>({l, r}) + shift, x);
            tMax.add(Range<ll>({l, r}) + shift, x);
            tMin.add(Range<ll>({l, r}) + shift, x);
        } else if (opt == 2) {
            // get min
            assert(*min_element(fakeSegtree + l, fakeSegtree + r + 1) == tMin.query(Range<ll>({l, r}) + shift));
        } else if (opt == 3) {
            // get max
            assert(*max_element(fakeSegtree + l, fakeSegtree + r + 1) == tMax.query(Range<ll>({l, r}) + shift));
        } else if (opt == 4) {
            // get sum
            function<bool(MinSegtree<ll>*)> isMin = [&](MinSegtree<ll>* node) {
                return node->query() == tMin.query();
            };
            assert(min_element(begin(fakeSegtree), end(fakeSegtree)) - begin(fakeSegtree) == tMin.findFirst(isMin)->l() + shift);
        } else if (opt == 5) {
            // get sum
            // function<bool(Segtree*)> isMax = [&](Segtree* node) {
            //     return node.query() == tMax.query();
            // };
            // assert(max_element(fakeSegtree + l, fakeSegtree + r + 1) - fakeSegtree == tMax.findFirst(isMax)->l());
        } else if (opt == 6) {
            int expected = 0;
            for (int i = l; i <= r; ++i) expected += fakeSegtree[i];
            assert(expected == tSum.query(Range<ll>({l, r}) + shift));
        }
    }
}

signed main() {
    srand(time(NULL));
    //testSegtree_manual();
    testSegtree_auto();
}
