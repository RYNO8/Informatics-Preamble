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

int main() {
    cin.tie(0); ios::sync_with_stdio(0);
    //ifstream cin{ "in.txt" };
    cin >> N >> A >> B >> C;
    for (ll i = 1; i <= N; ++i) cin >> arr[i];
    for (ll i = 1; i <= N; ++i) prefix[i] = prefix[i - 1] + arr[i];
    
    CHT cht;
    dp[0] = 0;
    for (ll i = 0; i <= N; ++i) {
        if (i == 0) dp[i] = 0;
        else dp[i] = cht.getMinima(prefix[i]) + A * prefix[i] * prefix[i] + B * prefix[i] + C; // Ax^2 + Bx;

        CHTLine l = {
            .m = -2 * A * prefix[i],
            .b = A * prefix[i] * prefix[i] - B * prefix[i] + dp[i]
        };
        cht.addLine(l);
    }
    cout << dp[N] << "\n";
#ifdef DEBUG
    cout << cht;
#endif
}