struct MinSt {
    ll l, r, m, minVal = 0;
    ll lazyVal = 0, lazyFlag = 0;
    MinSt* lChild = nullptr, * rChild = nullptr;
    MinSt(ll _l = 0, ll _r = MAXN - 1) {
        l = _l; r = _r;
        m = (l + r) / 2;
    }

    void doPush(MinSt* from, MinSt* to) {
        if (from->lazyFlag == 2) { // absolute, always override
            to->minVal = from->lazyVal;
            to->lazyVal = from->lazyVal;
            to->lazyFlag = 2;
        }
        else if (from->lazyFlag == 1) { // relative, add to absolute
            to->minVal += from->lazyVal;
            to->lazyVal += from->lazyVal;
            if (to->lazyFlag == 0) to->lazyFlag = 1;
        }
    }

    void push() {
        if (lChild == nullptr) lChild = new MinSt(l, m);
        if (rChild == nullptr) rChild = new MinSt(m + 1, r);
        doPush(this, lChild);
        doPush(this, rChild);
        lazyVal = 0;
        lazyFlag = 0;
    }

    void add(ll tl, ll tr, ll x) {
        // range update: add
        if (l > r || tl > r || l > tr || tl > tr) return;
        else if (tl <= l && l <= r && r <= tr) {
            minVal += x;
            lazyVal += x;
            if (lazyFlag == 0) lazyFlag = 1;
            return;
        }

        push();
        lChild->add(tl, min(m, tr), x);
        rChild->add(max(tl, m + 1), tr, x);
        minVal = min(lChild->minVal, rChild->minVal);
    }

    void set(ll tl, ll tr, ll x) {
        // range update: set
        if (l > r || tl > r || l > tr || tl > tr) return;
        else if (tl <= l && l <= r && r <= tr) {
            minVal = x;
            lazyVal = x;
            lazyFlag = 2;
            return;
        }

        push();
        lChild->set(tl, min(m, tr), x);
        rChild->set(max(tl, m + 1), tr, x);
        minVal = min(lChild->minVal, rChild->minVal);
    }

    ll getMin(ll tl, ll tr) {
        // range query: min
        if (l > r || tl > r || l > tr || tl > tr) return numeric_limits<ll>::max();
        else if (tl <= l && l <= r && r <= tr) return minVal;

        push();
        return min(lChild->getMin(tl, min(m, tr)), rChild->getMin(max(tl, m + 1), tr));
    }
};

struct MaxSt {
    ll l, r, m, maxVal = 0;
    ll lazyVal = 0, lazyFlag = 0;
    MaxSt* lChild = nullptr, * rChild = nullptr;
    MaxSt(ll _l = 0, ll _r = MAXN - 1) {
        l = _l; r = _r;
        m = (l + r) / 2;
    }

    void doPush(MaxSt* from, MaxSt* to) {
        if (from->lazyFlag == 2) { // absolute, always override
            to->maxVal = from->lazyVal;
            to->lazyVal = from->lazyVal;
            to->lazyFlag = 2;
        }
        else if (from->lazyFlag == 1) { // relative, add to absolute
            to->maxVal += from->lazyVal;
            to->lazyVal += from->lazyVal;
            if (to->lazyFlag == 0) to->lazyFlag = 1;
        }
    }

    void push() {
        if (lChild == nullptr) lChild = new MaxSt(l, m);
        if (rChild == nullptr) rChild = new MaxSt(m + 1, r);
        doPush(this, lChild);
        doPush(this, rChild);
        lazyVal = 0;
        lazyFlag = 0;
    }

    void add(ll tl, ll tr, ll x) {
        // range update: add
        if (l > r || tl > r || l > tr || tl > tr) return;
        else if (tl <= l && l <= r && r <= tr) {
            maxVal += x;
            lazyVal += x;
            if (lazyFlag == 0) lazyFlag = 1;
            return;
        }

        push();
        lChild->add(tl, min(m, tr), x);
        rChild->add(max(tl, m + 1), tr, x);
        maxVal = max(lChild->maxVal, rChild->maxVal);
    }

    void set(ll tl, ll tr, ll x) {
        // range update: set
        if (l > r || tl > r || l > tr || tl > tr) return;
        else if (tl <= l && l <= r && r <= tr) {
            maxVal = x;
            lazyVal = x;
            lazyFlag = 2;
            return;
        }

        push();
        lChild->set(tl, min(m, tr), x);
        rChild->set(max(tl, m + 1), tr, x);
        maxVal = max(lChild->maxVal, rChild->maxVal);
    }

    ll getMax(ll tl, ll tr) {
        // range query: max
        if (l > r || tl > r || l > tr || tl > tr) return 0;
        else if (tl <= l && l <= r && r <= tr) return maxVal;

        push();
        return max(lChild->getMax(tl, min(m, tr)), rChild->getMax(max(tl, m + 1), tr));
    }
};
struct SumSt {
    ll l, r, m, sumVal = 0;
    ll lazyVal = 0, lazyFlag = 0;
    SumSt* lChild = nullptr, * rChild = nullptr;
    SumSt(ll _l = 0, ll _r = MAXN - 1) {
        l = _l; r = _r;
        m = (l + r) / 2;
    }
    void doPush(SumSt* from, SumSt* to) {
        if (from->lazyFlag == 2) { // absolute, always override
            to->sumVal = from->lazyVal * (from->r - from->l + 1) / 2;
            to->lazyVal = from->lazyVal;
            to->lazyFlag = 2;
        }
        else if (from->lazyFlag == 1) { // relative, add to absolute
            to->sumVal += from->lazyVal * (from->r - from->l + 1) / 2;
            to->lazyVal += from->lazyVal;
            if (to->lazyFlag == 0) to->lazyFlag = 1;
        }
    }

    void push() {
        if (l == r) throw 20; // is leaf
        if (lChild == nullptr) lChild = new SumSt(l, m);
        if (rChild == nullptr) rChild = new SumSt(m + 1, r);
        doPush(this, lChild);
        doPush(this, rChild);
        lazyFlag = 0;
        lazyVal = 0;
    }

    void add(ll tl, ll tr, ll x) {
        // range update: add
        if (tl > r || l > tr) return;
        else if (tl <= l && r <= tr) {
            sumVal += x * (r - l + 1);
            lazyVal += x;
            if (lazyFlag == 0) lazyFlag = 1;
            return;
        }

        push();
        lChild->add(tl, min(m, tr), x);
        rChild->add(max(tl, m + 1), tr, x);
        sumVal = lChild->sumVal + rChild->sumVal;
    }

    void set(ll tl, ll tr, ll x) {
        // range update: set
        if (tl > r || l > tr) return;
        else if (tl <= l && l <= r && r <= tr) {
            sumVal = x * (r - l + 1);
            lazyVal = x;
            lazyFlag = 2;
            return;
        }

        push();
        lChild->set(tl, min(m, tr), x);
        rChild->set(max(tl, m + 1), tr, x);
        sumVal = lChild->sumVal + rChild->sumVal;
    }

    ll getSum(ll tl, ll tr) {
        // range query: max
        if (l > r || tl > r || l > tr || tl > tr) return 0LL;
        else if (tl <= l && l <= r && r <= tr) return sumVal;

        push();
        return lChild->getSum(tl, min(tr, m)) + rChild->getSum(max(tl, m + 1), tr);
    }
};

// Lazy range tree
// TODO: check
// TODO: persistency
class Segtree {
    /************************************************
     *                INITIALISATION                *
     ************************************************/

private:
    ll L, R;
    MinSt* minTree;
    MaxSt* maxTree;
    SumSt* totalTree;

    void checkRange(ll tl, ll tr) {
        // instead, maybe narrow the range until it is valid?
        assert(L <= tl && tl <= tr && tr <= R && "Invalid range");
    }

public:
    Segtree(ll _l = 0, ll _r = MAXN - 1) {
        assert(0 <= _l && _l <= _r && "Invalid range");
        L = _l;
        R = _r;
        minTree = new MinSt(L, R);
        maxTree = new MaxSt(L, R);
        totalTree = new SumSt(L, R);
    }

    /************************************************
     *                    DISPLAY                   *
     ************************************************/
     // [ O(1) ]
    void print(ostream& out = cout, bool newLine = true) {
        for (ll i = L; i <= min(L + 10, R); ++i) cout << get(i) << ' ';
        cout << '\n';

        if (newLine) cout << '\n';
    }

    /************************************************
     *                     UPDATES                  *
     ************************************************/
     // [ O(log N) ]
    void set(ll tl, ll tr, ll x) {
        checkRange(tl, tr);
        minTree->set(tl, tr, x);
        maxTree->set(tl, tr, x);
        totalTree->set(tl, tr, x);
    }
    void set(ll i, ll x) {
        set(i, i, x);
    }

    // [ O(log N) ]
    void clear() {
        set(L, R, 0);
    }

    // [ O(log N) ]
    void add(ll tl, ll tr, ll x) {
        checkRange(tl, tr);
        minTree->add(tl, tr, x);
        maxTree->add(tl, tr, x);
        totalTree->add(tl, tr, x);
    }
    void add(ll i, ll x) {
        add(i, i, x);
    }

    /************************************************
     *                    QUERIES                   *
     ************************************************/

    // [ O(log N) ]
    ll get(ll i) {
        checkRange(i, i);
        return getMin(i, i);
    }
    ll operator[](ll i) {
        return get(i);
    }

    // [ O(log N) ]
    ll getMin(ll tl, ll tr) {
        checkRange(tl, tr);
        return minTree->getMin(tl, tr);
    }
    ll getMin(ll i) {
        return getMin(i, i);
    }

    // [ O(log N) ]
    ll getMax(ll tl, ll tr) {
        checkRange(tl, tr);
        return maxTree->getMax(tl, tr);
    }
    ll getMax(ll i) {
        return getMax(i, i);
    }

    // [ O(log N) ]
    ll getTotal(ll tl, ll tr) {
        checkRange(tl, tr);
        return totalTree->getSum(tl, tr);
    }
    ll getTotal(ll i) {
        return getTotal(i, i);
    }

    /************************************************
     *                   ITERABLES                  *
     ************************************************/
     // TODO: pointers & std functionality
     // TODO: lower bound, upper bound, binary search (tree walk)
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
ostream& operator<<(ostream& out, Segtree tree) {
    tree.print(out);
    return out;
}