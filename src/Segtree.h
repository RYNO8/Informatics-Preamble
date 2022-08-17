#pragma once

// `T` should preferably be a real number
template<typename T> class Segtree {
    /************************************************
     *                INITIALISATION                *
     ************************************************/

private:
    ll l, r, m;
    T minVal, maxVal, sumVal;
    T lazyVal = T(0);
    short lazyMode = 0; // needs to hold {0, 1, 2}, is there a better way to do this?
    Segtree* lChild = nullptr, * rChild = nullptr;

    // O(1) push lazy values to children, for lazy propagation
    void push() {
        if (l == r) return;

        if (lChild == nullptr) lChild = new Segtree(l, m);
        if (rChild == nullptr) rChild = new Segtree(m + 1, r);

        if (lazyMode == 2) { // absolute, always override
            lChild->minVal = lazyVal;
            rChild->minVal = lazyVal;
            lChild->maxVal = lazyVal;
            rChild->maxVal = lazyVal;
            lChild->sumVal = lazyVal * T(m - l + 1);
            rChild->sumVal = lazyVal * T(r - (m + 1) + 1);

            lChild->lazyVal = lazyVal;
            rChild->lazyVal = lazyVal;
            lChild->lazyMode = 2;
            rChild->lazyMode = 2;
        }
        else if (lazyMode == 1) { // relative, add lChild absolute
            lChild->minVal += lazyVal;
            rChild->minVal += lazyVal;
            lChild->maxVal += lazyVal;
            rChild->maxVal += lazyVal;
            lChild->sumVal += lazyVal * T(m - l + 1);
            rChild->sumVal += lazyVal * T(r - (m + 1) + 1);

            lChild->lazyVal += lazyVal;
            rChild->lazyVal += lazyVal;
            if (lChild->lazyMode == 0) lChild->lazyMode = 1;
            if (rChild->lazyMode == 0) rChild->lazyMode = 1;
        }

        lazyVal = T(0);
        lazyMode = 0;
    }

    // O(1) update values from values of children
    void pull() {
        minVal = min(lChild->minVal, rChild->minVal);
        maxVal = max(lChild->maxVal, rChild->maxVal);
        sumVal = lChild->sumVal + rChild->sumVal;
    }

public:
    // O(1) Initialises an empty segtree containing the first `n` elements, which are all 0 
    Segtree(ll n) : l(0), r(n - 1) {
        assert(0 <= l && l <= r);
        m = (l + r) / 2;
        minVal = maxVal = sumVal = 0;
    }

    // O(1) Initialises an empty segtree containing elements from `l` to `r` inclusive, which are all 0 
    Segtree(ll _l, ll _r, T initVal = 0) : l(_l), r(_r) {
        assert(0 <= l && l <= r);
        m = (l + r) / 2;
        minVal = maxVal = sumVal = lazyVal = initVal;
        lazyMode = 2;
    }

    // O(1) Initialises segtree from existing segtree
    Segtree(Segtree<T>& t) {
        *this = t;
    }

    // O(N log N) Initialises segtree from existing segtree, in a the given range from `l` to `r` inclusive
    Segtree(Segtree<T> t, ll _l, ll _r) : l(_l), r(_r) {
        assert(0 <= l && l <= r);
        m = (l + r) / 2;
        for (ll i = _l; i <= _r; ++i) set(i, t.getVal(i));
    }

    // O(1) Initialise segtree given an existing left and right child
    Segtree(Segtree<T> _lChild, Segtree<T> _rChild) : l(_lChild.l), r(_rChild.r), lChild(&_lChild), rChild(&_rChild) {
        assert(0 <= l && l <= r);
        assert(lChild->r + 1 == rChild->l);
        m = (l + r) / 2;
        minVal = min(lChild->minVal, rChild->minVal);
        maxVal = max(lChild->maxVal, rChild->maxVal);
        sumVal = lChild->sumVal + rChild->sumVal;
    }

    inline bool covers(ll i) const {
        return l <= i && i <= r;
    }
    inline bool covers(ll tl, ll tr) const {
        return l <= tl && tl <= tr && tr <= r;
    }
    inline bool coveredBy(ll tl, ll tr) const {
        return tl <= l && r <= tr;
    }
    inline bool touching(ll tl, ll tr) const {
        return tl <= r && l <= tr && tl <= tr;
    }

    /************************************************
     *                    DISPLAY                   *
     ************************************************/

     // O(1) print the first 20 elements
    void print(ostream& out = cout, bool newLine = false) {
        for (ll i = l; i <= min(l + 20, r); ++i) out << getVal(i) << ' ';
        if (newLine) out << '\n';
    }

    // O(1)
    void printDebug(ostream& out = cout, bool newLine = false) {
        out << '[' << l << ".." << r << "] min=" << minVal << ", max=" << maxVal << ", sum=" << sumVal;
        if (newLine) out << '\n';
    }

    /************************************************
     *                     UPDATES                  *
     ************************************************/

     // O(log N) Range update: add
    Segtree* add(ll tl, ll tr, T x, bool persistent = false) {
        if (tl > r || l > tr) return this;

        push();
        Segtree* newNode = persistent ? new Segtree(*this) : this;
        if (coveredBy(tl, tr)) {
            newNode->minVal += x;
            newNode->maxVal += x;
            newNode->sumVal += x * T(r - l + 1);
            newNode->lazyVal += x;
            if (newNode->lazyMode == 0) newNode->lazyMode = 1;
        }
        else if (touching(tl, tr)) {
            newNode->lChild = lChild->add(tl, tr, x, persistent);
            newNode->rChild = rChild->add(tl, tr, x, persistent);
            newNode->pull();
        }
        return newNode;
    }
    Segtree* add(ll i, T x) {
        return add(i, i, x);
    }

    // O(log N) Range update: set
    Segtree* set(ll tl, ll tr, T x, bool persistent = false) {
        if (tl > r || l > tr) return this;

        push();
        Segtree* newNode = persistent ? new Segtree(*this) : this;
        if (coveredBy(tl, tr)) {
            newNode->minVal = x;
            newNode->maxVal = x;
            newNode->sumVal = x * T(r - l + 1);
            newNode->lazyVal = x;
            newNode->lazyMode = 2;
        }
        else if (touching(tl, tr)) {
            newNode->lChild = lChild->set(tl, tr, x, persistent);
            newNode->rChild = rChild->set(tl, tr, x, persistent);
            newNode->pull();
        }
        return newNode;
    }
    Segtree* set(ll i, T x) {
        return set(i, i, x);
    }

    // O(log N) Range update: set all values to T(0)
    void clear() {
        set(l, r, T(0));
    }

    /************************************************
     *                    QUERIES                   *
     ************************************************/

     // O(log N) Range query: mininum
    T getMin(ll tl, ll tr) {
        if (!touching(tl, tr)) return numeric_limits<T>::max();
        if (coveredBy(tl, tr)) return minVal;

        push();
        return min(lChild->getMin(tl, tr), rChild->getMin(tl, tr));
    }
    // O(1) Range query: min of all entries
    T getMin() {
        return getMin(l, r);
    }

    // O(log N) Range query: maximum
    T getMax(ll tl, ll tr) {
        if (!touching(tl, tr)) return numeric_limits<T>::min();
        if (coveredBy(tl, tr)) return maxVal;

        push();
        return max(lChild->getMax(tl, tr), rChild->getMax(tl, tr));
    }
    // O(1) Range query: max of all entries
    T getMax() {
        return getMax(l, r);
    }

    // O(log N) Range query: sum
    T getSum(ll tl, ll tr) {
        if (!touching(tl, tr)) return T(0); // additive identity
        else if (coveredBy(tl, tr)) return sumVal;

        push();
        return lChild->getSum(tl, tr) + rChild->getSum(tl, tr); // hoping theres no overflow
    }

    // O(1) Range query: sum of all entries
    T getSum() {
        return getSum(l, r);
    }
    

    ll getMinIndexFirst(ll tl, ll tr, T target) {
        if (!touching(tl, tr)) return -1;
        if (coveredBy(tl, tr) && target == minVal) {
            Segtree* node = this;
            while (node->l != node->r) {
                node->push();
                node = node->lChild->minVal == target ? node->lChild : node->rChild;
            }
            return node->l;
        }

        push(); // this is actually unnecessary, since `getMin` has already been called
        ll lIndex = lChild->getMinIndexFirst(tl, tr);
        if (lIndex != -1) return lIndex;
        return rChild->getMinIndexFirst(tl, tr);
    }
    // O(log N) Returns the index corresponding to the first minimum value in the range
    ll getMinIndexFirst(ll tl, ll tr) {
        T target = getMin(tl, tr);
        return getMinIndexFirst(tl, tr, target);
    }
    // O(log N) Returns the index corresponding to the first global minimum value
    ll getMinIndexFirst() {
        return getMinIndexFirst(l, r);
    }


    ll getMinIndexLast(ll tl, ll tr, T target) {
        if (!touching(tl, tr)) return -1;
        if (coveredBy(tl, tr) && target == minVal) {
            Segtree* node = this;
            while (node->l != node->r) {
                node->push();
                node = node->rChild->minVal == target ? node->rChild : node->lChild;
            }
            return node->l;
        }

        push(); // this is actually unnecessary, since `getMin` has already been called
        ll rIndex = rChild->getMinIndexLast(tl, tr);
        if (rIndex != -1) return rIndex;
        return lChild->getMinIndexLast(tl, tr);
    }
    // O(log N) Returns the index corresponding to the last minimum value in the range
    ll getMinIndexLast(ll tl, ll tr) {
        T target = getMin(tl, tr);
        return getMinIndexLast(tl, tr, target);
    }
    // O(log N) Returns the index corresponding to the last global minimum value
    ll getMinIndexLast() {
        return getMinIndexLast(l, r);
    }


    ll getMaxIndexFirst(ll tl, ll tr, T target) {
        if (!touching(tl, tr)) return -1;
        if (coveredBy(tl, tr) && target == maxVal) {
            Segtree* node = this;
            while (node->l != node->r) {
                node->push();
                node = node->lChild->maxVal == target ? node->lChild : node->rChild;
            }
            return node->l;
        }

        push(); // this is actually unnecessary, since `getMin` has already been called
        ll lIndex = lChild->getMaxIndexFirst(tl, tr);
        if (lIndex != -1) return lIndex;
        return rChild->getMaxIndexFirst(tl, tr);
    }
    // O(log N) Returns the index corresponding to the first maximum value in the range
    ll getMaxIndexFirst(ll tl, ll tr) {
        T target = getMax(tl, tr);
        return getMaxIndexFirst(tl, tr, target);
    }
    // O(log N) Returns the index corresponding to the first global maximum value
    ll getMaxIndexFirst() {
        return getMaxIndexFirst(l, r);
    }


    ll getMaxIndexLast(ll tl, ll tr, T target) {
        if (!touching(tl, tr)) return -1;
        if (coveredBy(tl, tr) && target == maxVal) {
            Segtree* node = this;
            while (node->l != node->r) {
                node->push();
                node = node->rChild->maxVal == target ? node->rChild : node->lChild;
            }
            return node->l;
        }

        push(); // this is actually unnecessary, since `getMin` has already been called
        ll rIndex = rChild->getMaxIndexLast(tl, tr);
        if (rIndex != -1) return rIndex;
        return lChild->getMaxIndexLast(tl, tr);
    }
    // O(log N) Returns the index corresponding to the first maximum value in the range
    ll getMaxIndexLast(ll tl, ll tr) {
        T target = getMax(tl, tr);
        return getMaxIndexLast(tl, tr, target);
    }
    // O(log N) Returns the index corresponding to the first global maximum value
    ll getMaxIndexLast() {
        return getMaxIndexLast(l, r);
    }

    // O(N log N) get the indexes of all elements with a minimum value
    vector<ll> getMinIndexes() {
        vector<ll> output;
        getMinIndexes(output);
        return output;
    }
    vector<ll> getMinIndexes(vector<ll>& output) {
        if (l == r) output.push_back(l);
        else {
            push();
            if (lChild->minVal == minVal) lChild->getMinIndexes(output);
            if (rChild->minVal == minVal) rChild->getMinIndexes(output);
        }
        return output;
    }

    // O(N log N) get the indexes of all elements with a maximum value
    vector<ll> getMaxIndexes() {
        vector<ll> output;
        getMaxIndexes(output);
        return output;
    }
    vector<ll> getMaxIndexes(vector<ll>& output) {
        if (l == r) output.push_back(l);
        else {
            push();
            if (lChild->maxVal == maxVal) lChild->getMaxIndexes(output);
            if (rChild->maxVal == maxVal) rChild->getMaxIndexes(output);
        }
        return output;
    }

    // O(log N) get value at index `i`
    T getVal(ll i) {
        assert(covers(i) && "Invalid index");
        push();
        if (l == r) {
            assert(minVal == maxVal && maxVal == sumVal && "AHH implementation broken");
            return minVal;
        }
        else if (i <= m) return lChild->getVal(i);
        else return rChild->getVal(i);
    }
    T operator[](ll i) {
        return getVal(i);
    }

};


/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Segtree<T> tree) {
    tree.print(out);
    return out;
}


