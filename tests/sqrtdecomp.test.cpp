// solution to Omo
// https://orac.amt.edu.au/cgi-bin/train/problem.pl?problemid=1102
/*

SAMPLE INPUT
22 6
87 79 87 95 83 85 67 72 95 68 79 71 69 95 71 79 68 95 72 67 85 83
4 22
9 9
10 17
10 11
9 17
10 18

SAMPLE OUTPUT
YES
YES
NO
NO
YES
YES

*/
#include "../src/Constants.h"
#include "../src/SqrtDecomp.h"
#include "../src/Ranges.h"
using namespace DS;
using namespace std;

void testSqrtDecomp() {
    int N, Q, bad = 0, arr[100005], curr[100005];
    vector<pair<int, int>> queries;

    cin >> N >> Q;
    for (int i = 0; i < N; ++i) cin >> arr[i];
    for (int a, b, i = 0; i < Q; ++i) {
        cin >> a >> b;
        queries.push_back({ a - 1, b - 1 });
    }

    coordCompressTransform(arr, arr + N);

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

    for (bool ans : chunkQueries<int, bool>(queries.begin(), queries.end(), add, rem, answer)) cout << (ans ? "YES" : "NO") << '\n';
}

signed main() {
    testSqrtDecomp();
}