#ifndef RANGES_H
#define RANGES_H
#include <algorithm>
#include <set>
#include <map>
#include <limits>
#include <type_traits>
#include <vector>
#include <iostream>
#include <assert.h>

namespace DS {
    // retrieve the underlying value given an iterator
    // e.g. `value_type<std::vector<int>::begin>` is `int`
    template<typename Iterator>
    using value_type = typename std::iterator_traits<Iterator>::value_type;

    // O(N log N)
    // @returns A mapping between original coordinates and compressed coordinates
    template<typename T, std::enable_if_t<std::is_integral_v<value_type<T>>, bool> = true>
    std::map<value_type<T>, value_type<T>> coordCompressMap(T begin, T end, int startI = 0) {
        std::vector<value_type<T>> unique(begin, end);
        sort(unique.begin(), unique.end());
        auto last = std::unique(unique.begin(), unique.end());
        unique.erase(last, unique.end());

        std::map<value_type<T>, value_type<T>> compressed;

        int i = startI;
        for (auto it = unique.begin(); it != unique.end(); ++it, ++i) {
            compressed[*it] = i;
        }
        return compressed;
    }

    // Compresses coordinates in place
    template<typename T, std::enable_if_t<std::is_integral_v<value_type<T>>, bool> = true>
    void coordCompressTransform(T begin, T end, int startI = 0){
        std::map<value_type<T>, value_type<T>> m = coordCompressMap(begin, end, startI);
        for (T it = begin; it != end; it++) *it = m[*it];
    }

    // ahhh its probably bad to inherit from std types, but i really want the comparison operators ;-;
    template<
        typename T,
        std::enable_if_t<std::is_integral_v<T> || std::is_floating_point_v<T>, bool> = true
    >
    struct Range : public std::pair<T, T> {
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

        using std::pair<T, T>::pair;
        using std::pair<T, T>::first;
        using std::pair<T, T>::second;

        // inclusive inclusive
        Range(T lVal, T rVal): std::pair<T, T>(lVal, rVal) {
            assert(l() <= r());
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        friend std::ostream& operator<<(std::ostream &out, const Range &range) {
            out << '[' << range.l() << ".." << range.r() << ']';
            return out;
        }

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        // O(1)
        // @returns Inclusive left endpoint of range
        T l() const {
            return first;
        }

        // O(1)
        // @returns Inclusive right endpoint of range
        T r() const {
            return second;
        }

        // O(1)
        // @returns Midpoint of range (rounding down)
        constexpr T midpoint() const {
            // NOTE: using the implementation so it still rounds down when
            return l() + (r() - l()) / 2;
        }

        // O(1)
        // @returns Inclusive inclusive length of range
        // @note Has support for integral and non integral types T
        constexpr T length() const {
            return r() - l() + (std::numeric_limits<T>::is_integer ? T(1) : T(0));
        }

        // O(1)
        // @returns
        bool isEmpty() const {
            return l() == r();
        }

        /************************************************
         *                  OPERATIONS                  *
         ************************************************/

        // O(1)
        // @returns
        bool covers(T i) const {
            return l() <= i && i <= r();
        }
        
        // O(1)
        // @returns
        bool covers(const Range<T> &o) const {
            return l() <= o.l() && o.r() <= r();
        }
        

        // O(1)
        // @returns
        bool coveredBy(const Range<T> &o) const {
            return o.l() <= l() && r() <= o.r();
        }

        // O(1)
        // @returns
        bool touches(const Range<T> &o) const  {
            return l() == o.r() || o.l() == r();
        }

        // O(1)
        // @returns
        bool overlaps(const Range<T> &o) const  {
            return l() <= o.r() && o.l() <= r();
        }

        bool canMerge(const Range<T> &o) const {
            return l() <= o.r() + 1 && o.l() <= r() + 1;
        }

        // what to do if range is invalid?
        // O(1) intersect
        // @returns
        Range<T> operator&(const Range<T> &o) const {
            return { std::max(l(), o.l()), std::min(r(), o.r()) };
        }

        Range<T> operator|(const Range<T> &o) const {
            assert(canMerge(o));
            return { std::min(l(), o.l()), std::max(r(), o.r()) };
        }

        // should be defined by comparison on pairs

        // bool operator<(const Range<T> o) const {
        //     if (l == o.l) return r < o.r;
        //     return l < o.l;
        // }

        // bool operator>(const Range<T> o) const {
        //     if (l == o.l) return r > o.r;
        //     return l > o.l;
        // }

        // bool operator<=(const Range<T> o) const {
        //     if (l == o.l) return r <= o.r;
        //     return l <= o.l;
        // }

        // bool operator>=(const Range<T> o) const {
        //     if (l == o.l) return r >= o.r;
        //     return l >= o.l;
        // }

        /************************************************
         *               TRANSFORMATIONS                *
         ************************************************/

        Range<T> operator+(T shift) {
            return { l() + shift, r() + shift };
        }

        Range<T> operator-(T shift) {
            return { l() - shift, r() - shift };
        }

        Range<T> operator*(T scale) {
            return { l() * scale, r() * scale };
        }

        std::pair<Range<T>, Range<T>> split(T val) {

        }
    };

    template<typename T, std::enable_if_t<std::is_integral_v<T>, bool> = true>
    class Ranges: public std::set<Range<T>> {
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

    public:
        // inherit constructors and all constructor overloads
        using std::set<Range<T>>::set;
        using std::set<Range<T>>::begin;
        using std::set<Range<T>>::end;
        using std::set<Range<T>>::rbegin;
        using std::set<Range<T>>::rend;
        using std::set<Range<T>>::lower_bound;
        using std::set<Range<T>>::upper_bound;

        // void _cleanup() {
        //     // merge ranges when
        // }

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        const T length() const {
            T total = T(0);
            for (Range<T> r : *this) {
                total += r.length();
            }
            return total;
        }

        /************************************************
         *                  OPERATIONS                  *
         ************************************************/

        // O(log N)
        // @returns
        bool covers(int i) const {
            auto it = lower_bound(Range<T>(i + 1, i + 1));
            if (it == begin()) return false;
            it--;
            return it->covers(i);
        }

        // O(log N)
        // @returns
        bool covers(const Range<T> &o) const {
            auto it = lower_bound(Range<T>(o.l() + 1, o.l() + 1));
            if (it == begin()) return false;
            it--;
            return it->covers(o);
        }

        // O(N)
        // @returns
        bool coversAll(const Ranges<T> &o) const {
            return (*this & o) == o;
        }

        // O(N log N)
        // intersect
        // @returns
        Ranges<T> operator&(const Ranges<T> &o) const {
            Ranges<T> out;
            auto it2 = o.begin();
            for (auto it1 = begin(); it1 != end(); ++it1) {
                while (it2 != o.end() && it2->r() < it1->l()) ++it2;
                while (it2 != o.end() && it1->overlaps(*it2)) {
                    out.insert(*it1 & *it2);
                    ++it2;
                }
                --it2;
            }
            return out;
        }

        // O(N log N)
        // union
        Ranges<T> operator|(const Ranges<T> &o) {
            std::vector<Range<T>> out;
            for (auto it1 = begin(), it2 = o.begin(); it1 != end() || it2 != o.end(); ) {
                auto it = (it2 == o.end() || (it1 != end() && it1->l() <= it2->l())) ? it1++ : it2++;
                if (!out.empty() && out.back().canMerge(*it)) {
                    out.back() = out.back() | *it;
                } else {
                    out.push_back(*it);
                }
            }
            return Ranges<T>(out.begin(), out.end());
        }

        // O(N log N)
        Ranges<T> operator~() const {
            Ranges<T> out;
            T prevEnd = std::numeric_limits<T>::min();
            for (Range<T> r : *this) {
                if (prevEnd <= r.l() - 1) out.insert({ prevEnd, r.l() - 1 });
                prevEnd = r.r() + 1;
            }
            if (rbegin()->r() + 1 <= std::numeric_limits<T>::max()) {
                out.insert({ rbegin()->r() + 1, std::numeric_limits<T>::max() });
            }
            return out;
        }

        bool isDisjoint(const Ranges<T> &o) const {
            return (*this & o).empty();
        }

        /************************************************
         *               TRANSFORMATIONS                *
         ************************************************/

        // O(N log N)
        // `ranges.split(on)` should be equivalent to `ranges & ~Ranges({on})`
        Ranges<T> split(Range<T> &on) {
            Ranges<T> out;
            for (Range<T> r : *this) {
                if (r.covers(on)) {
                    std::pair<Range<T>, Range<T>> partition = r.split(on);
                    if (!partition.first.isEmpty()) out.insert(partition.first);
                    if (!partition.second.isEmpty()) out.insert(partition.second);
                } else {
                    out.insert(r);
                }
            }
            return out;
        }
    };
};

#endif
