#include "../src/Constants.h"
#include "../src/SqrtDecomp.h"
using namespace std;

void testSqrtDecomp() {
	int N, Q, bad = 0, arr[100005], curr[100005];
	vector<pair<int, int>> queries;

	cin >> N >> Q;
	for (int i = 1; i <= N; ++i) cin >> arr[i];
	for (int a, b, i = 0; i < Q; ++i) {
		cin >> a >> b;
		queries.push_back({ a, b });
	}

	coordCompressTransform(arr + 1, arr + N + 1);

	function<void(int)> add = [&](int i) {
		if (curr[arr[i]] % 2 == 0) ++bad;
		else --bad;
		curr[arr[i]]++;
	};

	function<void(int)> rem = [&](int i) {
		if (curr[arr[i]] % 2 == 0) ++bad;
		else --bad;
		curr[arr[i]]--;
	};

	function<bool(int, int)> answer = [&](int l, int r) {
		return bad == ((r - l + 1) % 2 == 1);
	};

	for (bool ans : chunkQueries(queries, add, rem, answer)) cout << (ans ? "YES" : "NO") << '\n';
}

signed main() {
	testSqrtDecomp();
}