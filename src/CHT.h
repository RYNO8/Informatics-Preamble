#ifndef CHT_H
#define CHT_H
#include "Constants.h"
#include "Util.h"

namespace DS {
// Representation of an infinite line using intersection gradient form
// TODO: bad, see CHT in icpc template
struct CHTLine {
    ll m, b;

    // O(1)
    // Displays the line
    // @param `out` The string representation of the graph is piped to this output stream
    // @param `newLine` Indicates whether to end with a trailing `\\n`
    friend std::ostream &operator<<(std::ostream &out, const CHTLine &line) {
        out << "y = ";
        if (line.m != 0) {
            out << line.m << 'x';
            if (line.b >= 0) out << '+';
        }
        out << line.b;
        return out;
    }

    // O(1)
    // @returns the x value of the intersection
    double intersect(CHTLine a) const {
        return (double)(b - a.b) / (a.m - m);
    }

    // O(1)
    // @returns Evalute the y value when x = x
    // @note May overflow
    ll operator()(ll x) const {
        return m * x + b;
    }
};

// Convex Hull Trick
// "infinte" range and domain, but added lines need to have non-decreasing gradient

class CHT {
public:

    friend std::ostream &operator<<(std::ostream &out, const CHT &c) {
        return out << c.L;
    }

    // amortised O(1)
    // Include a unique line into the convex hull, where `m` is non-decreasing
    void addLine(CHTLine l) {
        while (L.size() >= 2 && L.back().intersect(*++L.rbegin()) >= L.back().intersect(l)) {
            L.pop_back();
        }
        // while (L.size() >= 2) {
        //     CHTLine b = L.back();
        //     CHTLine a = *++L.rbegin();
        //     if ((a.b - b.b) * (l.m - b.m) >= (b.b - l.b) * (b.m - a.m)) {
        //         L.pop_back();
        //     } else break;
        // }
        while (!L.empty() && L.back().m == l.m && L.back().b <= l.b) {
            L.pop_back();
        }
        if (L.empty() || l.m > L.back().m || l.b > L.back().b) L.push_back(l);
    }

    // O(log n)
    // @returns Minimum y value when x=`x`
    ll getMaxima(ll x) const {
        if (L.empty()) return LLONG_MIN;
        int l = -1, r = L.size() - 1;
        while (l + 1 < r) {
            int m = (l + r) / 2;
            // if (L[m].b - L[m+1].b >= x * (L[m+1].m - L[m].m)) r = m;
            if (L[m].intersect(L[m + 1]) >= x) r = m;
            else l = m;
        }
        return L[r](x);
    }

private:
    std::vector<CHTLine> L;
};

// Li Chao tree
// Lines can be inserted and queried in any order, but range is limited to `MAXN`
// Initially contains the line 0x+0
class LiChaoTree {
private:
    CHTLine tree[2 * MAXN];

public:
    // O(1)
    // Initialisation
    LiChaoTree(CHTLine initLine) {
        std::fill(std::begin(tree), std::end(tree), initLine);
    }

    // O(log MAXN)
    // Inserts a line
    void addLine(CHTLine line, int node = 1, int l = 0, int r = MAXN - 1) {
        int m = (l + r) / 2;
        bool lGood = line(l) < tree[node](l);
        bool mGood = line(m) < tree[node](m);
        if (mGood) std::swap(tree[node], line);

        if (l == r) return;
        else if (lGood != mGood) addLine(line, node * 2, l, m);
        else addLine(line, node * 2 + 1, m + 1, r);
    }

    // O(log MAXN)
    // @returns Minimum y value when x=`x`
    ll getMinima(int x, int node = 1, int l = 0, int r = MAXN - 1) {
        assert(l <= x && x <= r && "x out of range");
        int m = (l + r) / 2;
        ll ans = tree[node](x);
        if (l == r) {}
        else if (x <= m) ans = std::min(ans, getMinima(x, node * 2, l, m));
        else ans = std::min(ans, getMinima(x, node * 2 + 1, m + 1, r));
        return ans;
    }
};
}; // namespace DS

#endif
