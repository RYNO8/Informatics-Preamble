#pragma once
#include "Constants.h"

namespace DS {
    // Representation of an infinite line using intersection gradient form
    struct CHTLine {
        ll m, b;

        // O(1)
        // Displays the line
        // @param `out` The string representation of the graph is piped to this output stream
        // @param `newLine` Indicates whether to end with a trailing `\\n`
        void print(std::ostream& out = std::cout, bool newLine = true) const {
            out << "y = " << m << " x + " << b;
            if (newLine) out << '\n';
        }

        // O(1)
        // @returns the x value of the intersection
        double intersect(CHTLine a) const {
            return (double)(b - a.b) / (a.m - m);
        }

        // O(1)
        // @returns Evalute the y value when x = x
        // @note May overflow
        ll eval(ll x) const {
            return m * x + b;
        }
    };

    std::ostream& operator<<(std::ostream& out, CHTLine line) {
        line.print(out);
        return out;
    }

    // Convex Hull Trick
    // "infinte" range and domain, but added lines need to have non-decreasing gradient
    class CHT : public std::vector<CHTLine> {
public:
        // amortised O(1)
        // Include a unique line into the convex hull, where `m` is non-decreasing
        void addLine(CHTLine l) {
            while (size() >= 2 && back().intersect(*++rbegin()) >= back().intersect(l)) {
                pop_back();
            }
            while (!empty() && back().m == l.m && back().b >= l.b) {
                pop_back();
            }
            push_back(l);
        }

        // O(log n)
        // @returns Minimum y value when x=`x`
        ll getMinima(ll x) const {
            if (empty()) return std::numeric_limits<ll>::max();
            int l = -1, r = size() - 1;
            while (l + 1 < r) {
                int m = (l + r) / 2;
                if ((*this)[m].intersect((*this)[m + 1]) >= x) r = m;
                else l = m;
            }
            return (*this)[r].eval(x);
        }
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
        LiChaoTree() {

        }

        // O(log MAXN)
        // Inserts a line
        void addLine(CHTLine line) {
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
            // implemented iteratively
            int node = 1, l = 0, r = MAXN - 1;
            ll ans = LLONG_MIN;
            assert(l <= x && x <= r && "x out of range");

            while (true) {
                ans = std::max(ans, tree[node].eval(x));
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
};