#pragma once
#include "Constants.h"
#include "Ranges.h"

namespace DS {
    // `T` should preferably be a real number
    // however i dont want to use "requires std::integral<T> || std::floating_point<T>"
    // because segtree should support ModInt

    // I should be left and right identity of Combine
    // Combine should have `T operator()(T, T)` which is associative property
    // CombineAgg should have `operator()(T x, unsigned ll n)` which is defined as adding x to itself n times (n > 0)
    // TODO: if not provided use fast exponentiation
    // Null should be the left (and right?) identity of Update
    // Update should have `T operator()(T, T)`
    // TODO: what properties required?
    template<typename T, T I, class Combine, class CombineAgg, T Null, class Update>
    class Segtree : public Range<ll> {
        /************************************************
         *                INITIALISATION                *
         ************************************************/

    private:
        enum LazyMode {
            // no lazy needs to be performed here
            CLEAN = 0,

            // relative, corresponds to `data := Update(data, v)` operations (i.e. augmenting operations)
            AUGMENT = 1,

            // absolute (always override), corresponds to `data := v` operations
            SET = 2
        };

        T val = I;
        T lazy = Null;
        LazyMode lazyMode = CLEAN;
        Segtree *lChild = nullptr, *rChild = nullptr;

        bool isLeaf() const {
            return isEmpty();
        }

        // O(1) push lazy values to children, for lazy propagation
        void push() {
            if (isLeaf()) return;

            if (lChild == nullptr) lChild = new Segtree(l(), midpoint());
            if (rChild == nullptr) rChild = new Segtree(midpoint() + 1, r());

            if (lazyMode == SET) { // absolute, always override
                lChild->val = CombineAgg()(lazy, (ull)lChild->length());
                rChild->val = CombineAgg()(lazy, (ull)rChild->length());

                lChild->lazy = lazy;
                rChild->lazy = lazy;
                lChild->lazyMode = SET;
                rChild->lazyMode = SET;
            }
            else if (lazyMode == AUGMENT) { // relative, add lChild absolute
                lChild->val = Update()(lChild->val, CombineAgg()(lazy, (ull)lChild->length()));
                rChild->val = Update()(rChild->val, CombineAgg()(lazy, (ull)rChild->length()));

                lChild->lazy = Update()(lChild->lazy, lazy);
                rChild->lazy = Update()(rChild->lazy, lazy);
                if (lChild->lazyMode == CLEAN) lChild->lazyMode = AUGMENT;
                if (rChild->lazyMode == CLEAN) rChild->lazyMode = AUGMENT;
            }

            lazy = Null;
            lazyMode = CLEAN;
        }

        // O(1) update values from values of children
        void pull() {
            val = Combine()(lChild->val, rChild->val);
        }

    public:
        // O(1) Initialises an empty segtree containing elements from `l` to `r` inclusive, which are all I
        Segtree(ll _l, ll _r, T initVal = I) : Range(_l, _r) {
            val = CombineAgg()(initVal, (ull)length());
            lazy = initVal;
            lazyMode = SET;
        }

        // O(N log N) Initialises segtree from existing segtree, in a the given range from `l` to `r` inclusive
        Segtree(Segtree<T, I, Combine, CombineAgg, Null, Update> t, ll _l, ll _r) : Range(_l, _r) {
            for (ll i = l(); i <= r(); ++i) set(i, t.query(i));
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        // O(1) print the first 20 elements
        void print(std::ostream& out = std::cout, bool newLine = false) {
            for (ll i = l(); i <= std::min(l() + 20, r()); ++i) out << query(i) << ' ';
            if (newLine) out << '\n';
        }

        // O(1)
        void printDebug(std::ostream& out = std::cout, bool newLine = false) const {
            Range::print();
            out << " = " << val << '\n';
            if (newLine) out << '\n';
        }

        /************************************************
         *                     UPDATES                  *
         ************************************************/

        // O(log N) Range update: add
        // performs T[l] = Update(T[l], x), ..., T[r] = Update(T[r], x),
        // @returns the root of the updated segtree
        Segtree* add(Range<ll> queryRange, T x, bool persistent = false) {
            if (!overlaps(queryRange)) return this;

            
            Segtree* newNode = persistent ? new Segtree(*this) : this;
            if (coveredBy(queryRange)) {
                newNode->val = Update()(newNode->val, CombineAgg()(x , (ull)length()));
                newNode->lazy = Update()(newNode->lazy, x);
                if (newNode->lazyMode == CLEAN) newNode->lazyMode = AUGMENT;
            }
            else {
                push();
                newNode->lChild = lChild->add(queryRange, x, persistent);
                newNode->rChild = rChild->add(queryRange, x, persistent);
                newNode->pull();
            }
            return newNode;
        }
        Segtree* add(ll i, T x) {
            return add(i, i, x);
        }

        // O(log N) Range update: set
        // performs T[l] = x, ..., T[r] = x
        // @returns the root of the updated segtree
        Segtree* set(Range<ll> queryRange, T x, bool persistent = false) {
            if (!overlaps(queryRange)) return this;

            
            Segtree* newNode = persistent ? new Segtree(*this) : this;
            if (coveredBy(queryRange)) {
                newNode->val = CombineAgg()(x, (ull)length());
                newNode->lazy = x;
                newNode->lazyMode = SET;
            } else {
                push();
                newNode->lChild = lChild->set(queryRange, x, persistent);
                newNode->rChild = rChild->set(queryRange, x, persistent);
                newNode->pull();
            }
            return newNode;
        }
        Segtree* set(ll i, T x) {
            return set({i, i}, x, false);
        }
        Segtree* set(T x) {
            return set({l, r}, x, false);
        }

        // O(log N) Range update: set all values to I
        void clear() {
            set({l, r}, I);
        }

        /************************************************
         *                    QUERIES                   *
         ************************************************/

        // O(log N) Range query
        // @returns Combine(T[l], Combine(..., Combine(T[r-1], T[r])...))
        // @returns I if the segtree range has no overlap with the query range
        T query(Range<ll> queryRange) {
            if (!overlaps(queryRange)) return I;
            else if (coveredBy(queryRange)) return val;

            push();
            return Combine()(lChild->query(queryRange), rChild->query(queryRange)); // hoping theres no overflow
        }
        T operator[](Range<ll> queryRange) {
            return query(queryRange);
        }

        // O(1) Range query: sum of all entries
        T query() {
            return query(Range(l(), r()));
            // or
            // return val;
        }
        // O(log N) get value at index `i`
        T query(ll i) {
            assert(covers(i) && "Invalid index");
            push();
            if (l == r) return val;
            else if (i <= midpoint()) return lChild->query(i);
            else return rChild->query(i);
        }
        T operator[](ll i) {
            return query(i);
        }

        // treewalk down the filtered tree, finding the leftmost lowest node
        Segtree* findFirst(std::function<bool(Segtree*)> isTarget) {
            if (!isTarget(this)) return nullptr;
            Segtree* node = this;
            while (!node->isLeaf()) {
                node->push();
                if (isTarget(node->lChild)) node = node->lChild;
                else if (isTarget(node->rChild)) node = node->rChild;
                else return node;
            }
            return node;
        }

        Segtree* findLast(std::function<bool(Segtree*)> isTarget) {
            if (!isTarget(this)) return nullptr;
            Segtree* node = this;
            while (!node->isLeaf()) {
                node->push();
                if (isTarget(node->rChild)) node = node->rChild;
                else if (isTarget(node->lChild)) node = node->lChild;
                else return node;
            }
            return node;
        }

        Segtree* extendRight(int initVal) {
            Segtree *root = new Segtree(l(), r() + length());
            root->lazy = Null;
            root->lazyMode = CLEAN;
            root->lChild = this;
            root->rChild = new Segtree(l() + length(), r() + length(), initVal);
            root->pull();
            return root;
        }

        Segtree* extendLeft(int initVal) {
            Segtree *root = new Segtree(l(), r() + length());
            root->lazy = Null;
            root->lazyMode = CLEAN;
            root->lChild = new Segtree(l() + length(), r() + length(), initVal);
            root->rChild = this;
            root->pull();
            return root;
        }
    };

    template<typename T> class MaxOp {
public:
        constexpr T operator()(const T &lhs, const T &rhs) const {
            return std::max(lhs, rhs);
        }
    };

    template<typename T> class MinOp {
public:
        constexpr T operator()(const T &lhs, const T &rhs) const {
            return std::min(lhs, rhs);
        }
    };

    template<typename T> class FirstOp {
public:
        constexpr T operator()(const T &lhs, const ll &rhs) const {
            return lhs;
        }
    };

    template<typename T>
    using SumSegtree = Segtree<T, 0, std::plus<T>, std::multiplies<T>, 0, std::plus<T>>;

    template<typename T>
    using MaxSegtree = Segtree<T, std::numeric_limits<T>::min(), MaxOp<T>, FirstOp<T>, 0, std::plus<T>>;

    template<typename T>
    using MinSegtree = Segtree<T, std::numeric_limits<T>::max(), MinOp<T>, FirstOp<T>, 0, std::plus<T>>;

    /************************************************
     *                    DISPLAY                   *
     ************************************************/
    template<typename T, T I, class Combine, class CombineAgg, T Null, class Update>
    std::ostream& operator<<(std::ostream& out, Segtree<T, I, Combine, CombineAgg, Null, Update> tree) {
        tree.print(out);
        return out;
    }
};
