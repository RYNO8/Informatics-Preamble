#ifndef SEGTREE_H
#define SEGTREE_H
#include "Constants.h"
#include "Ranges.h"
#include <functional>
#include <algorithm>
#include <type_traits>

namespace DS {
    // `T` should preferably be a real number
    // because segtree should support ModInt

    // I want to support compile time specification of whether Segtree is persistent, but
    // error: non-type template parameters of class type only available with ‘-std=c++2a’ or ‘-std=gnu++2a’
    // template<std::enable_if_t<!is_persistent::value, bool> = true>
    // thankfully c++ metaprogramming has a true_type
    // template<typename is_persistent_ = is_persistent, std::enable_if_t<!is_persistent_::value, bool> = true>
    // template<typename is_persistent_ = is_persistent, std::enable_if_t<is_persistent_::value, bool> = true>
    // oop turns out I'm not gonna do constexpr stuff
    
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
        typename Combine,
        typename CombineAgg,
        T Null,
        typename Update,
        typename is_persistent,
        std::enable_if_t<is_my_integral_v<T>, bool> = true
    >
    class Segtree : public Range<ll> {
        using MySegtree = Segtree<T, I, Combine, CombineAgg, Null, Update, is_persistent>;

        /************************************************
         *                INITIALISATION                *
         ************************************************/

    private:
        enum LazyMode {
            // no lazy needs to be performed here
            CLEAN,

            // relative, corresponds to `data := Update(data, v)` operations (i.e. augmenting operations)
            AUGMENT,

            // absolute (always override), corresponds to `data := v` operations
            SET
        };

        T val = I;
        T lazy = Null;
        LazyMode lazyMode = CLEAN;
        MySegtree *lChild = nullptr, *rChild = nullptr, *parent = nullptr;

        bool isLeaf() const {
            return isEmpty();
        }

        MySegtree* copy() const {
            return new MySegtree(*this);
        }

        // O(1) push lazy values to children, for lazy propagation
        void push() {
            if (isLeaf()) return;

            // if (lChild == nullptr) lChild = new MySegtree(l(), midpoint());
            // else if (is_persistent::value) lChild = lChild->copy();
            // if (rChild == nullptr) rChild = new MySegtree(midpoint() + 1, r());
            // else if (is_persistent::value) rChild = rChild->copy();

            if (lChild == nullptr) {
                lChild = new MySegtree(l(), midpoint());
                rChild = new MySegtree(midpoint() + 1, r());
            } else if (is_persistent::value) {
                lChild = lChild->copy();
                rChild = rChild->copy();
            }

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
        Segtree(ll _l, ll _r, MySegtree t) : Range(_l, _r) {
            for (ll i = l(); i <= r(); ++i) set(i, t.query(i));
        }

        // O(N log N)
        // Initialise from brace enclosed initializer list
        Segtree(T* begin, T* end) : Range(0, std::distance(begin, end) - 1) {
            size_t i = 0;
            for (auto it = begin; it != end; it++, ++i) set(i, *it);
        }

        // O(N log N)
        // Initialise from brace enclosed initializer list
        Segtree(typename std::vector<T>::iterator begin, typename std::vector<T>::iterator end) : Range(0, std::distance(begin, end) - 1) {
            size_t i = 0;
            for (auto it = begin; it != end; it++, ++i) set(i, *it);
        }

        // O(N log N)
        // Initialise from brace enclosed initializer list
        template<size_t N>
        Segtree(typename std::array<T, N>::iterator begin, typename std::array<T, N>::iterator end) : Range(0, std::distance(begin, end) - 1) {
            size_t i = 0;
            for (auto it = begin; it != end; it++, ++i) set(i, *it);
        }

        // O(N)
        // Create the entire tree (such that is it no longer lazy propagation)
        // Cleans all nodes (no more lazy values)
        // This could help with constant factor optimisations, since initailising nodes is expensive
        // @note Modifies tree in place
        void pushAll() {
            if (!isLeaf()) {
                push();
                lChild->pushAll();
                rChild->pushAll();
            }
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        // O(N)
        // @note: not const segtree because queries do pushes
        friend std::ostream& operator<<(std::ostream &out, Segtree segtree) {
            out << '[';
            for (ll i = segtree.l(); i <= segtree.r(); ++i) {
                out << ' ' << segtree.query(i);
            }
            out << " ]";
            return out;
        }

        // O(displayRange.size())
        // @note: not const segtree because queries do pushes
        // @TODO this a butchered solution until I do range views
        std::ostream& display(std::ostream &out, Range<ll> displayRange) {
            out << '[';
            for (ll i = displayRange.l(); i <= r(); ++i) {
                out << ' ' << query(i);
            }
            out << " ]";
            return out;
        }

        // O(1)
        // @note first part looks like display of Range<ll>, but this code has been copied
        friend std::ostream& operator<<(std::ostream &out, const Segtree* segtree) {
            out << '[' << segtree->l() << ".." << segtree->r() << "] = " << segtree->val;
            return out;
        }

        

        /************************************************
         *                     UPDATES                  *
         ************************************************/

        // // O(log N) Range update: add
        // // performs T[l] = Update(T[l], x), ..., T[r] = Update(T[r], x),
        // // @note Modifies segtree in place
        void add(Range<ll> queryRange, T x) {
            if (!overlaps(queryRange)) return;
            else if (coveredBy(queryRange)) {
                val = Update()(val, CombineAgg()(x , (ull)length()));
                lazy = Update()(lazy, x);
                if (lazyMode == CLEAN) lazyMode = AUGMENT;
            }
            else {
                push();
                lChild->add(queryRange, x);
                rChild->add(queryRange, x);
                pull();
            }
        }

        // O(log N) Point update: add
        // performs T[i] = Update(T[i], x)
        void add(ll i, T x) {
            add({i, i}, x);
        }

        // O(log N) Range update: set
        // performs T[l] = x, ..., T[r] = x
        // @note Modifies segtree in place
        void set(const Range<ll> &queryRange, T x) {
            if (!overlaps(queryRange)) return;
            else if (coveredBy(queryRange)) {
                val = CombineAgg()(x, (ull)length());
                lazy = x;
                lazyMode = SET;
            } else {
                push();
                lChild->set(queryRange, x);
                rChild->set(queryRange, x);
                pull();
            }
        }

        // O(log N) Point update: set
        // performs T[i] = x
        void set(ll i, T x) {
            set({i, i}, x);
        }

        // O(log N) Range update: set all values to I
        void clear() {
            set({l(), r()}, I);
        }

        /************************************************
         *                    QUERIES                   *
         ************************************************/

        // O(log N) Range query: Combine
        // @returns Combine(T[l], Combine(..., Combine(T[r-1], T[r])))
        // @returns I if the segtree range has no overlap with the query range
        T query(const Range<ll> &queryRange) {
            if (!overlaps(queryRange)) return I;
            else if (coveredBy(queryRange)) return val;

            push();
            return Combine()(lChild->query(queryRange), rChild->query(queryRange));
        }

        // TODO: pandas view kinda thing?
        T operator[](const Range<ll> queryRange) {
            return query(queryRange);
        }

        // O(1) Range query: sum of all entries
        T query() const {
            return val;
        }
        // O(log N) get value at index `i`
        T query(ll i) {
            assert(covers(i) && "Invalid index");
            if (isLeaf()) return val;
            push();
            if (i <= midpoint()) return lChild->query(i);
            else return rChild->query(i);
        }

        T operator[](ll i) {
            return query(i);
        }

    public:
        // treewalk down the filtered tree, finding the leftmost lowest node
        MySegtree* findFirst(std::function<bool(MySegtree*)> isTarget) {
            if (!isTarget(this)) return nullptr;
            MySegtree* node = this;
            while (!node->isLeaf()) {
                node->push();
                if (isTarget(node->lChild)) node = node->lChild;
                else if (isTarget(node->rChild)) node = node->rChild;
                else return node;
            }
            return node;
        }

        MySegtree* findLast(std::function<bool(MySegtree*)> isTarget) {
            if (!isTarget(this)) return nullptr;
            MySegtree* node = this;
            while (!node->isLeaf()) {
                node->push();
                if (isTarget(node->rChild)) node = node->rChild;
                else if (isTarget(node->lChild)) node = node->lChild;
                else return node;
            }
            return node;
        }

        MySegtree* extendRight(int initVal) {
            MySegtree* root = new MySegtree(l(), r() + length());
            root->lazy = Null;
            root->lazyMode = CLEAN;
            root->lChild = this;
            root->rChild = new MySegtree(l() + length(), r() + length(), initVal);
            root->pull();
            return root;
        }

        MySegtree* extendLeft(int initVal) {
            MySegtree* root = new MySegtree(l(), r() + length());
            root->lazy = Null;
            root->lazyMode = CLEAN;
            root->lChild = new MySegtree(l() + length(), r() + length(), initVal);
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
        MySegtree* findLeaf(ll i) {
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
            using value_type        = MySegtree const;
            using pointer           = MySegtree const*;
            using reference         = MySegtree const&;
            using internal_type     = MySegtree*;

        protected:
            ll pos;
            internal_type m_ptr;

        public:
            iterator(MySegtree *segtree, ll pos_) : pos(pos_) {
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
    using SumSegtree = Segtree<T, 0, std::plus<T>, std::multiplies<T>, 0, std::plus<T>, std::false_type>;

    template<typename T>
    using SumSegtreePersistent = Segtree<T, 0, std::plus<T>, std::multiplies<T>, 0, std::plus<T>, std::true_type>;


    template<typename T>
    using MaxSegtree = Segtree<T, std::numeric_limits<T>::min(), MaxOp<T>, FirstOp<T>, 0, std::plus<T>, std::false_type>;

    template<typename T>
    using MaxSegtreePersistent = Segtree<T, std::numeric_limits<T>::min(), MaxOp<T>, FirstOp<T>, 0, std::plus<T>, std::true_type>;


    template<typename T>
    using MinSegtree = Segtree<T, std::numeric_limits<T>::max(), MinOp<T>, FirstOp<T>, 0, std::plus<T>, std::false_type>;

    template<typename T>
    using MinSegtreePersistent = Segtree<T, std::numeric_limits<T>::max(), MinOp<T>, FirstOp<T>, 0, std::plus<T>, std::true_type>;

};

#endif