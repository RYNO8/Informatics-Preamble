#pragma once
#include "Constants.h"

// Representation of an infinite line using intersection gradient form
struct Line_ {
    ll m, b;

    // O(1)
    // Displays the line
    // @param `out` The string representation of the graph is piped to this output stream
    // @param `newLine` Indicates whether to end with a trailing `\\n`
    void print(ostream& out = cout, bool newLine = true) const {
        out << m << " x + " << b;
        if (newLine) out << '\n';
    }

    // O(1)
    // @returns the x value of the intersection
    double intersect(Line_ a) const {
        return (double)(b - a.b) / (a.m - m);
    }

    // O(1)
    // @returns Evalute the y value when x = x
    // @note May overflow
    ll eval(ll x) const {
        return m * x + b;
    }
};

ostream& operator<<(ostream& out, Line_ line) {
    line.print(out);
    return out;
}

// Convex Hull Trick
// "infinte" range and domain, but added lines need to have non-decreasing gradient
struct CHT {
    vector<Line_> cht;

    // O(1)
    // Initialisation
    CHT() {
        cht.clear();
    }

    // O(N)
    // Displays the convex hull, showing the lines in increasing order
    // @param `out` The string representation of the graph is piped to this output stream
    // @param `newLine` Indicates whether to end with a trailing `\\n`
    void print(ostream& out = cout, bool newLine = true) const {
        out << cht;
        if (newLine) out << '\n';
    }

    // amortised O(1)
    // Include a unique line into the convex hull, where `m` is non-decreasing
    void addLine(Line_ l) {
        while (cht.size() >= 2 && cht[cht.size() - 1].intersect(cht[cht.size() - 2]) >= cht[cht.size() - 1].intersect(l)) {
            cht.pop_back();
        }
        while (!cht.empty() && cht.back().m == l.m && cht.back().b >= l.b) {
            cht.pop_back();
        }
        cht.push_back(l);
    }

    // O(log n)
    // @returns Minimum y value when x=`x`
    ll getLine(ll x) const {
        if (cht.empty()) return numeric_limits<ll>::max();
        int l = -1, r = cht.size() - 1;
        while (l + 1 < r) {
            int m = (l + r) / 2;
            if (cht[m].intersect(cht[m + 1]) >= x) r = m;
            else l = m;
        }
        return cht[r].eval(x);
    }
};

ostream& operator<<(ostream& out, CHT cht) {
    cht.print(out);
    return out;
}

// Li Chao tree
// Lines can be inserted and queried in any order, but range is limited to `MAXN`
// Initially contains the line 0x+0
struct LiChaoTree {
    Line_ tree[2 * MAXN];

    // O(1)
    // Initialisation
    LiChaoTree() {

    }

    // O(log MAXN)
    // Inserts a line
    void addLine(Line_ line) {
        int node = 1, l = 0, r = MAXN - 1;
        do {
            int m = (l + r) / 2;
            bool lGood = line.eval(l) >= tree[node].eval(l);
            bool mGood = line.eval(m) >= tree[node].eval(m);
            if (mGood) tree[node] = line;

            if (lGood != mGood) {
                node = node * 2;
                r = m;
            }
            else {
                node = node * 2 + 1;
                l = m + 1;
            }
        } while (l != r);
    }

    // O(log MAXN)
    // @returns Maximum y value when x=`x`
    ll getLine(int x) {
        int node = 1, l = 0, r = MAXN - 1;
        ll ans = LLONG_MIN;
        assert(l <= x && x <= r && "x out of range");

        while (true) {
            ans = max(ans, tree[node].eval(x));
            int m = (l + r) / 2;
            
            if (l == r) return ans;
            else if (x < m) {
                node = node * 2;
                r = m;
            }
            else {
                node = node * 2 + 1;
                l = m + 1;
            }
        }
    }
};