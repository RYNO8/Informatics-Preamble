struct Line {
    ll m, b;
    double intersect(Line a) {
        return (double)(b - a.b) / (a.m - m);
    }
    ll getVal(ll x) {
        return m * x + b;
    }
};

struct CHT {
    vector<Line> cht;
    CHT() {

    }

    // m should be non-decreasing
    void addLine(ll m, ll b) {
        Line l = { m, b };
        while (cht.size() >= 2 && cht[cht.size() - 1].intersect(cht[cht.size() - 2]) >= cht[cht.size() - 1].intersect(l)) {
            cht.pop_back();
        }
        while (!cht.empty() && cht.back().m == m && cht.back().b >= b) {
            cht.pop_back();
        }
        cht.push_back(l);
    }

    // returns minimum line intersection
    ll getLine(ll x) {
        if (cht.empty()) return numeric_limits<ll>::max();
        int l = -1, r = cht.size() - 1;
        while (l + 1 < r) {
            int m = (l + r) / 2;
            if (cht[m].intersect(cht[m + 1]) >= x) r = m;
            else l = m;
        }
        return cht[r].getVal(x);
    }
};