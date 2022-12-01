#include <assert.h>
#include <limits>
#include <set>
#include <vector>
#pragma once

namespace DS {
    // retrieve the underlying value given an iterator
    // e.g. `value_type<std::vector<int>::begin>` is `int`
    template<typename Iterator>
    using value_type = typename std::iterator_traits<Iterator>::value_type;

    // @returns A mapping between original coordinates and compressed coordinates
    template<typename T> std::map<value_type<T>, value_type<T>> coordCompressMap(T begin, T end, int startI = 0) {
        std::map<value_type<T>, value_type<T>> compressed;

        int i = startI;
        for (auto val = begin; val != end; ++val, ++i) {
            compressed[*val] = i;
        }
        return compressed;
    }
    // Compresses coordinates in place
    template<typename T> void coordCompressTransform(T begin, T end, int startI = 0){
        std::map<value_type<T>, value_type<T>> m = coordCompressMap(begin, end, startI);
        for (T i = begin; i != end; i++) *i = m[*i];
    }

    // TODO: test
    // ahhh its probably bad to inherit from std types, but i really want the comparison operators ;-;
    template<typename T> struct Range : public std::pair<T, T> {
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

        // inclusive inclusive
        Range(T lVal, T rVal): std::pair<T, T>(lVal, rVal) {
            assert(l() <= r());
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        std::ostream& print(std::ostream &out = std::cout, bool newLine=false) const {
            out << '[';
            out << l() << ".." << r();
            out << ']';
            if (newLine) out << '\n';
            return out;
        }

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        // O(1)
        // @returns
        T l() const {
            // TODO: why cant i just do `return first;`
            return this->first;
        }

        // O(1)
        // @returns
        T r() const {
            // TODO: why cant i just do `return second;`
            return this->second;
        }

                // O(1)
                // @returns
        constexpr T midpoint() const {
            // NOTE: using the implementation so it still rounds down when
            return l() + (r() - l()) / 2;
        }

        // O(1)
        // @returns
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
        bool covers(const Range<T> &o) const {
            return l() <= o.l() && o.r() <= r();
        }
        // O(1)
        // @returns
        bool covers(T i) {
            return l() <= i <= r();
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

        // what to do if range is invalid?
        // O(1)
        // @returns
        Range<T> intersect(const Range<T> &o) const {
            return {max(l(), o.l()), min(r(), o.r())};
        }

        // O(1)
        // @returns
        Range<T> operator&(const Range<T> &o) const {
            return intersect(o);
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
    
    /************************************************
	 *                    DISPLAY                   *
	 ************************************************/
	template<typename T> std::ostream& operator<<(std::ostream& out, const Range<T>& val) {
		val.print(out);
		return out;
	}

    // TODO: this is currently a stub
    // TODO: test
    template<typename T> class Ranges: public std::set<Range<T>> {
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

    private:

    public:
        // inherit constructors and all constructor overloads
        using std::set<Range<T>>::set;

        // NOTE: dont need print because c++ will call print on the underlying set

        void _cleanup() {
            // merge ranges when
        }

        /************************************************
         *                  OPERATIONS                  *
         ************************************************/

        // should call this covers?
        bool contains(const Ranges<T> &o) const {
            return true;
        }

        bool overlaps(const Ranges<T> &o) const  {
            return true;
        }

        Ranges<T> getIntersect(const Ranges<T> &o) const {
            return true;
        }

        Ranges<T> operator*(const Ranges<T> &o) const {
            return getIntersect(o);
        }

        Ranges<T> getUnion(const Ranges<T> &o) const {
            return true;
        }

        Ranges<T> operator+(const Ranges<T> &o) const {
            return getUnion(o);
        }

        Ranges<T> getComplement(T rangeStart = std::numeric_limits<T>::min(), T rangeEnd = std::numeric_limits<T>::max()) const {

        }
        Ranges<T> operator~() const {
            return getComplement();
        }

        bool isDisjoint() const {
            return true;
        }

        Ranges<T> flatten() const {

        }

        /************************************************
         *               TRANSFORMATIONS                *
         ************************************************/

        Range<T> split(T val) {

        }

        Range<T> splits(std::vector<T> val) {

        }
    };

    /************************************************
	 *                    DISPLAY                   *
	 ************************************************/
	// template<typename T> std::ostream& operator<<(std::ostream& out, const Ranges<T>& val) {
	// 	val.print(out);
	// 	return out;
	// }
};