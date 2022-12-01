#include "../src/Segtree.h"
using namespace std;
using namespace DS;

void testSegtree_manual() {
	ll N;
	cin >> N;

	Segtree<ll> t(N);
	cout << "Enter operation: set (S), add (A), min (m), max (M), total (T), quit (Q)\n\n";

	while (true) {
		char op;
		cin >> op;

		ll tl, tr, x;
		switch (op) {
		case 'S':
			cin >> tl >> tr >> x;
			t = *t.set(tl, tr, x);
			break;
		case 'A':
			cin >> tl >> tr >> x;
			t = *t.add(tl, tr, x);
			break;
		case 'm':
			cin >> tl >> tr;
			cout << t.getMin(tl, tr) << '\n';
			break;
		case 'M':
			cin >> tl >> tr;
			cout << t.getMax(tl, tr) << '\n';
			break;
		case 'T':
			cin >> tl >> tr;
			cout << t.getSum(tl, tr) << '\n';
			break;
		case 'Q':
			return;
		}

		cout << t << '\n';

		ll iMin = t.getMinIndexFirst();
		assert(t[iMin] == t.getMin() && t.getMin(0, iMin - 1) > t[iMin]);
		cout << string(int(iMin) * 2, ' ') << "^ min\n";

		ll iMax = t.getMaxIndexFirst();
		assert(t[iMax] == t.getMax() && t.getMax(0, iMax - 1) < t[iMax]);
		cout << string(int(iMax) * 2, ' ') << "^ max\n";
	}
}


void testSegtree_auto() {
	int fakeSegtree[10000] = {};
	Segtree<int> t(10000);

	for (int rep = 0; rep < 10000; ++rep) {
		int opt = rand() % 5;
		int l = rand() % 10000, r = rand() % 10000;
		if (l > r) swap(l, r);

		if (opt == 0) {
			// set
			int x = (rand() % 100) - 50;
			for (int i = l; i <= r; ++i) fakeSegtree[i] = x;
			t.set(l, r, x);
		} else if (opt == 1) {
			// add
			int x = (rand() % 100) - 50;
			for (int i = l; i <= r; ++i) fakeSegtree[i] += x;
			t.add(l, r, x);
		} else if (opt == 2) {
			// get min
		} else if (opt == 3) {
			// get max
		} else {
			// get sum
			int expected = 0;
			for (int i = l; i <= r; ++i) expected += fakeSegtree[i];
			assert(expected == t.getSum(l, r));
		}
	}
}
signed main() {
	srand(time(NULL));
	//testSegtree_manual();
	testSegtree_auto();
}