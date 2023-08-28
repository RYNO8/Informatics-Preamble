// solution to Commando
// https://dmoj.ca/problem/apio10p1
/*
SAMPLE INPUT
4
-1 10 -20
2 2 3 4

SAMPLE OUTPUT
9
*/

#include "../src/CHT.h"
#include "../src/Util.h"
using namespace std;
using namespace DS;

ll N, A, B, C;
ll arr[1000005], prefix[1000005], dp[1000005];

CHT commando(istream &in, ostream &out) {
    in >> N >> A >> B >> C;
    for (ll i = 1; i <= N; ++i)
        in >> arr[i];
    for (ll i = 1; i <= N; ++i)
        prefix[i] = prefix[i - 1] + arr[i];

    CHT cht;
    dp[0] = 0;
    for (ll i = 0; i <= N; ++i) {
        if (i == 0)
            dp[i] = 0;
        else
            dp[i] = cht.getMinima(prefix[i]) + A * prefix[i] * prefix[i] + B * prefix[i] +
                    C; // Ax^2 + Bx;

        CHTLine l = {
            .m = -2 * A * prefix[i], .b = A * prefix[i] * prefix[i] - B * prefix[i] + dp[i]};
        cht.addLine(l);
    }
    out << dp[N] << "\n";
    return cht;
}

int main() {
    stringstream in, out;
    in << R"""(
4
-1 10 -20
2 2 3 4
)""";
    CHT cht = commando(in, out);
    assert(repr(cht) == "[ y = 0 y = 8x-52 y = 14x-114 y = 22x-222 ]");
    string ans;
    out >> ans;
    assert(ans == "9");
}
