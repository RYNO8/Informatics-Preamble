#include "../src/Constants.h"
#include "../src/Segtree.h"
using namespace std;

void testSegtree() {
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
signed main() {
	cin.tie(0); ios::sync_with_stdio(0);
	testSegtree();
}