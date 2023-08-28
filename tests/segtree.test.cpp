#include "../src/Segtree.h"
#include "../src/Constants.h"
#include "../src/Ranges.h"
#include "../src/Util.h"
using namespace std;
using namespace DS;

// solution to Lazy Update
// https://orac2.info/problem/datalazyupdate/
/*

SAMPLE INPUT
3 8
A 0 1 3
T 0 2
S 1 2 4
T 0 2
M 0 0
M 0 1
A 1 1 1
T 0 2

SAMPLE OUTPUT
6
11
3
4
12

*/
void testSegtree_manual() {
    ll N, Q;
    cin >> N >> Q;

    SumSegtree<ll> *tSum = new SumSegtree<ll>(0, N - 1, 0);
    MaxSegtree<ll> *tMax = new MaxSegtree<ll>(0, N - 1, 0);
    MinSegtree<ll> *tMin = new MinSegtree<ll>(0, N - 1, 0);
    // cout << "Enter operation: set (S), add (A), min (m), max (M), total (T), quit (Q)\n\n";

    for (int q = 0; q < Q; ++q) {
        char op;
        cin >> op;

        ll tl, tr, x;
        switch (op) {
        case 'S':
            cin >> tl >> tr >> x;
            tSum->set({tl, tr}, x);
            tMax->set({tl, tr}, x);
            tMin->set({tl, tr}, x);
            break;
        case 'A':
            cin >> tl >> tr >> x;
            tSum->add({tl, tr}, x);
            tMax->add({tl, tr}, x);
            tMin->add({tl, tr}, x);
            break;
        case 'm':
            cin >> tl >> tr;
            cout << tMin->query({tl, tr}) << '\n';
            break;
        case 'M':
            cin >> tl >> tr;
            cout << tMax->query({tl, tr}) << '\n';
            break;
        case 'T':
            cin >> tl >> tr;
            cout << tSum->query({tl, tr}) << '\n';
            break;
        }

        // cout << *tSum << '\n';

        // function<bool(MinSegtree<ll>*)> isMin = [&](MinSegtree<ll>* node) {
        //     return node->query() == tMin->query();
        // };
        // function<bool(MaxSegtree<ll>*)> isMax = [&](MaxSegtree<ll>* node) {
        //     return node->query() == tMax->query();
        // };

        // ll iMin = tMin->findFirst(isMin)->l();
        // cout << string(int(iMin) * 2, ' ') << "^ min\n";

        // ll iMax = tMax->findFirst(isMax)->l();
        // cout << string(int(iMax) * 2, ' ') << "^ max\n";
    }
}

struct fakeSegtree {
    int N;
    vector<int> t;
    fakeSegtree(int _N, int initVal) : N(_N) { t = vector<int>(N, initVal); }

    void set(int l, int r, int x) {
        for (int i = l; i <= r; ++i)
            t[i] = x;
    }

    void add(int l, int r, int x) {
        for (int i = l; i <= r; ++i)
            t[i] += x;
    }

    int sum(int l, int r) { return accumulate(t.begin() + l, t.begin() + r + 1, 0); }

    int min(int l, int r) { return *min_element(t.begin() + l, t.begin() + r + 1); }

    int max(int l, int r) { return *max_element(t.begin() + l, t.begin() + r + 1); }

    int minIndex() { return min_element(t.begin(), t.end()) - t.begin(); }

    int maxIndex() { return max_element(t.begin(), t.end()) - t.begin(); }

    void extendRight(int x) {
        int s = t.size();
        for (int rep = 0; rep < s; ++rep) {
            t.push_back(x);
        }
    }

    size_t size() { return t.size(); }
};

void testSegtree_auto() {
    srand(time(NULL));
    for (int tI = 0; tI < 20; ++tI) {
        // cout << "REP: " << tI << "\n";

        int N       = 5;
        int initVal = 69;
        fakeSegtree t(N, initVal);
        ll shift             = -rand();

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
                function<bool(MinSegtree<ll> *)> isMin = [&](MinSegtree<ll> *node) {
                    return node->query() == tMin->query();
                };
                assert(t.minIndex() == tMin->findFirst(isMin)->l() - shift);

            } else if (opt == 6) {
                // get index of max
                function<bool(MaxSegtree<ll> *)> isMax = [&](MaxSegtree<ll> *node) {
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

void testSegtree_benchmark() {
    srand(SIXTYNINE);
    SumSegtree<ll> *t = new SumSegtree<ll>(0, 100005);
    t->pushAll();
    for (int q = 0; q < 100000; ++q) {
        int opt = rand() % 2;
        int l = rand() % 100000, r = rand() % 100000;

        if (opt == 0) {
            if (l > r) swap(l, r);
            t->set({l, r}, rand());
        } else {
            t->query({l, r});
        }
    }
}

void testSegtree_iterators() {
    ll a[]            = {1, 2, 3, 4};
    SumSegtree<ll> t_ = SumSegtree<ll>(begin(a), end(a));
    assert(repr(t_) == "[ 1 2 3 4 ]");
    assert(repr(&t_) == "[0..3] = 10");
    assert(repr(SumSegtree<ll>(1, 2, t_)) == "[ 2 3 ]");

    vector<ll> v{1, 2, 3, 4};
    SumSegtree<ll> t = SumSegtree<ll>(v.begin(), v.end());
    assert(repr(t) == "[ 1 2 3 4 ]");

    using It = SumSegtree<ll>::iterator;

    It it1   = t.begin();
    It it2   = ++t.begin();
    It it3   = ++ ++t.begin();
    It it4   = ++ ++ ++t.begin();
    It it5   = ++ ++ ++ ++t.begin();

    assert((*it1).query() == 1 && it1->query() == 1);
    assert((*it2).query() == 2 && it2->query() == 2);
    assert((*it3).query() == 3 && it3->query() == 3);
    assert((*it4).query() == 4 && it4->query() == 4);
    assert(it5 == t.end());

    assert(distance(t.begin(), t.end()) == 4);
    assert(distance(t.rbegin(), t.rend()) == 4);

    // satisfies DefaultConstructible
    // https://en.cppreference.com/w/cpp/named_req/DefaultConstructible
    // default initialised
    It temp1;
    // value initialised
    It temp2{};
    // temporary
    It temp3 = It();
    It temp4 = It{};

    temp1    = it1;
    assert(--temp1 == t.end() /* && ++temp1 == it1*/);
    temp2 = it2;
    assert(--temp2 == it1 && ++temp2 == it2);
    temp3 = it3;
    assert(--temp3 == it2 && ++temp3 == it3);
    temp4 = it4;
    assert(--temp4 == it3 && ++temp4 == it4);

    assert(++t.rbegin() == it3);
    assert(t.rbegin() == it4);
    assert(--++t.rbegin() == it4);
    assert(t.rend() == t.end());

    assert(t.begin() + 6969 == t.end());
    assert(t.begin() - 6969 == t.end());
    assert(++t.end() == t.end());
    assert(2 + it1 == it3);
    assert(t.end() - t.begin() == 4);
    assert(t.begin()[2] == *it3);
    assert(t.end() > t.begin());
    assert(t.begin() >= t.begin());
    assert(t.rend() > t.rbegin());
}

signed main() {
    // testSegtree_manual();
    // testSegtree_iterators();
    testSegtree_benchmark();
    // testSegtree_auto();
}
