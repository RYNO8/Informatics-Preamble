#pragma once
#include "Constants.h"
#include "Ranges.h"
#include <functional>
#include <algorithm>

namespace DS {
    // `T` should preferably be a real number
    // because segtree should support ModInt

    // I should be left and right identity of Combine
    // Combine should have `T operator()(T, T)` which is associative property
    // CombineAgg should have `operator()(T x, unsigned ll n)` which is defined as adding x to itself n times (n > 0)
    // TODO: if not provided use fast exponentiation
    // Null should be the left (and right?) identity of Update
    // Update should have `T operator()(T, T)`
    // TODO: what properties required?
    template<
        typename T,
        T I,
        class Combine,
        class CombineAgg,
        T Null,
        class Update,
        std::enable_if_t<is_my_integral_v<T>, bool> = true
    >
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
        Segtree *lChild = nullptr, *rChild = nullptr, *parent = nullptr;

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
        // O(1)
        // Initialises an empty segtree containing elements from `l` to `r` inclusive, which are all I
        Segtree(ll _l, ll _r, T initVal = I) : Range(_l, _r) {
            val = CombineAgg()(initVal, (ull)length());
            lazy = initVal;
            lazyMode = SET;
        }

        // O(N log N)
        // Initialises segtree from existing segtree, in a the given range from `l` to `r` inclusive
        Segtree(ll _l, ll _r, Segtree t) : Range(_l, _r) {
            for (ll i = l(); i <= r(); ++i) set(i, t.query(i));
        }

        // O(N log N)
        // Initialise from brace enclosed initializer list
        Segtree(ll _l, ll _r, std::vector<T> vals) : Range(_l, _r) {
            size_t i = 0;
            for (ll pos = l(); i < vals.size() && pos <= r(); ++i, ++pos) set(pos, vals[i]);
        }

        // O(N log N)
        // Initialise from brace enclosed initializer list
        Segtree(std::vector<T> vals) : Range(0, vals.size() - 1) {
            size_t i = 0;
            for (ll pos = l(); i < vals.size() && pos <= r(); ++i, ++pos) set(pos, vals[i]);
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        // O(N)
        // note: not const segtree because queries do pushes
        friend std::ostream& operator<<(std::ostream& out, Segtree segtree) {
            out << "[ ";
            for (ll i = segtree.l(); i <= segtree.r(); ++i) {
                out << segtree.query(i);
                if (i != segtree.r()) out << ' ';
            }
            out << " ]";
            return out;
        }

        // O(1)
        friend std::ostream& operator<<(std::ostream& out, const Segtree* segtree) {
            out << Range<ll>(segtree) << " = " << segtree.val << '\n';
            return out;
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
            return add({i, i}, x);
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
        T query() const {
            return val;
        }
        // O(log N) get value at index `i`
        T query(ll i) {
            assert(covers(i) && "Invalid index");
            if (l() == r()) return val;
            push();
            if (i <= midpoint()) return lChild->query(i);
            else return rChild->query(i);
        }
        T operator[](ll i) {
            return query(i);
        }

    public:
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
            Segtree* root = new Segtree(l(), r() + length());
            root->lazy = Null;
            root->lazyMode = CLEAN;
            root->lChild = this;
            root->rChild = new Segtree(l() + length(), r() + length(), initVal);
            root->pull();
            return root;
        }

        Segtree* extendLeft(int initVal) {
            Segtree* root = new Segtree(l(), r() + length());
            root->lazy = Null;
            root->lazyMode = CLEAN;
            root->lChild = new Segtree(l() + length(), r() + length(), initVal);
            root->rChild = this;
            root->pull();
            return root;
        }


        /************************************************
         *                   ITERATOR                   *
         ************************************************/

        // TODO: test
        // TODO: test with persistent segtree
    private:
        enum ItPos {
            TOO_FAR_LEFT = -1,
            GOOD = 0,
            TOO_FAR_RIGHT = 1,
        };

        // @note almost identity to query(ll)
        Segtree* findLeaf(ll i) {
            assert(covers(i) && "Invalid index");
            if (l() == r()) return this;
            push();
            if (i <= midpoint()) {
                lChild->parent = this;
                return lChild->findLeaf(i);
            } else {
                rChild->parent = this;
                return rChild->findLeaf(i);
            }
        }

    public:
        // to do interator incr we need to backtrack up the tree
        // I have a choice to:
        // 1. store a backtrack chain of length log N in the iterator (log N space per iterator)
        // 2. parent pointers for each node in the segtree (constant additional space per node)
        // 3. variation of 2 - lazily generate parent pointers when iterating through tree
        // i have decided 3 and hopefully i wont regret it

        struct iterator {
            using iterator_category = std::random_access_iterator_tag;
            using difference_type   = std::ptrdiff_t;
            using value_type        = Segtree const;
            using pointer           = Segtree const*;
            using reference         = Segtree const&;
            using internal_type     = Segtree*;

        protected:
            ll pos;
            internal_type m_ptr;

        public:
            iterator(Segtree *segtree, ll pos_) : pos(pos_) {
                if (segtree->covers(pos)) {
                    m_ptr = segtree->findLeaf(pos);
                } else {
                    m_ptr = nullptr;
                }
            }
            iterator() : m_ptr(nullptr) {}
            iterator(internal_type ptr) : m_ptr(ptr) {
                assert(m_ptr == nullptr || m_ptr->isLeaf());
            }

            // O(1)
            // defeference the interator
            // @returns reference to the leaf node corresponding to the index
            reference operator*() const { 
                assert(m_ptr != nullptr && "cannot dereference out of bound interator");
                return *m_ptr;
            }

            // O(1)
            // access member of underlying of iterator
            // @returns pointer to the leaf node corresponding to the index
            pointer operator->() const {
                assert(m_ptr != nullptr && "cannot dereference out of bound interator");
                return m_ptr;
            }

            // O(log N)
            // Increment the iterator `n` times, by walking up towarsd the root
            // then walking down a subtree to the target
            iterator& operator+=(difference_type n) {
                pos += n;
                if (m_ptr != nullptr) {
                    ll target = m_ptr->l() + n;
                    while (m_ptr != nullptr && !m_ptr->covers(target)) {
                        m_ptr = m_ptr->parent;
                        if (m_ptr != nullptr) m_ptr->pull();
                    }
                    if (m_ptr != nullptr) {
                        m_ptr = m_ptr->findLeaf(target);
                    }
                }
                return *this;
            }
            friend iterator operator+(const iterator it, difference_type n) { iterator temp = it; return temp += n; }
            friend iterator operator+(difference_type n, const iterator it) { iterator temp = it; return temp += n; }

            // O(log N)
            // Increment
            iterator& operator++() { return *this += 1; }

            // O(log N)
            // Postfix increment
            iterator operator++(int) { iterator tmp = *this; ++(*this); return tmp; }

            // O(log N)
            // Increment the iterator `n` times
            iterator& operator-=(difference_type n) { return *this += -n; }
            friend iterator operator-(const iterator it, difference_type n) { iterator temp = it; return temp -= n; }

            // O(1)
            // difference between 2 iterator positions
            friend difference_type operator-(const iterator it1, const iterator it2) {
                return it1.pos - it2.pos;
            }

            // O(log N)
            // Decrement
            iterator& operator--() { return *this -= 1; }

            // O(log N)
            // Postfix increment
            iterator operator--(int) { iterator tmp = *this; --(*this); return tmp; }

            // O(log N)
            // i[n] equivalent to *(i + n)
            reference operator[](difference_type n) const {
                return *(*this + n);
            }

            friend bool operator==(const iterator& a, const iterator& b) { return a.m_ptr == b.m_ptr; };
            friend bool operator!=(const iterator& a, const iterator& b) { return a.m_ptr != b.m_ptr; };
            // TODO: what if a and b are referencing different segtrees
            friend bool operator<(const iterator& a, const iterator& b) { return a.pos < b.pos; };
            friend bool operator>(const iterator& a, const iterator& b) { return a.pos > b.pos; };
            friend bool operator<=(const iterator& a, const iterator& b) { return a.pos <= b.pos; };
            friend bool operator>=(const iterator& a, const iterator& b) { return a.pos >= b.pos; };
        };

        struct reverse_iterator : public iterator {
            using iterator_category = typename std::iterator_traits<iterator>::iterator_category;
            using difference_type   = typename std::iterator_traits<iterator>::difference_type;
            using value_type        = typename std::iterator_traits<iterator>::value_type;
            using pointer           = typename std::iterator_traits<iterator>::pointer;
            using reference         = typename std::iterator_traits<iterator>::reference;
            using iterator::iterator;

            reverse_iterator& operator++() {
                iterator::operator--();
                return *this;
            }

            iterator& operator--() {
                iterator::operator++();
                return *this;
            }

            // // O(1)
            // difference between 2 iterator positions
            friend difference_type operator-(const reverse_iterator it1, const reverse_iterator it2) {
                return it2.pos - it1.pos;
            }

            // TODO: what if a and b are referencing different segtrees
            friend bool operator<(const reverse_iterator& a, const reverse_iterator& b) { return -a.pos < -b.pos; };
            friend bool operator>(const reverse_iterator& a, const reverse_iterator& b) { return -a.pos > -b.pos; };
            friend bool operator<=(const reverse_iterator& a, const reverse_iterator& b) { return -a.pos <= -b.pos; };
            friend bool operator>=(const reverse_iterator& a, const reverse_iterator& b) { return -a.pos >= -b.pos; };
        };

        // O(log N)
        iterator begin() {
            return iterator(this, l());
        }

        // O(1)
        iterator end() {
            return iterator(this, r() + 1);
        }

        // O(log N)
        reverse_iterator rbegin() {
            return reverse_iterator(this, r());
        }

        // O(1)
        reverse_iterator rend() {
            return reverse_iterator(this, l() - 1);
        }
    };

    /************************************************
     *                CUSTOM SEGTREES               *
     ************************************************/

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
};
