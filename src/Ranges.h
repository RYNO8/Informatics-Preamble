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
    template<typename T> class Range {
        /************************************************
         *                 INITIALISATION               *
         ************************************************/

    private:
        T lVal, rVal;
        bool lInclusive, rInclusive;

    public:
        Range(T lVal, T rVal, bool lInclusive=true, bool rInclusive=true): 
            lVal(lVal),
            rVal(rVal),
            lInclusive(lInclusive),
            rInclusive(rInclusive) {
            assert(lVal <= rVal);
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        std::ostream& print(std::ostream &out = std::cout, bool newLine=false) const {
            out << (lInclusive ? '[' : '(');
            out << lVal << ", " << rVal;
            out << (rInclusive ? ']' : ')');
            if (newLine) out << '\n';
            return out;
        }

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        const T& l() const {
            return &lVal;
        }

        const T& r() const {
            return &rVal;
        }

        T length() const {
            return (
                rVal - 
                lVal - 
                (rInclusive ? std::numeric_limits<T>::lowest() : 0) - 
                (lInclusive ? std::numeric_limits<T>::lowest() : 0)
            );
        }

        bool isNull() const {
            return length() == 0;
        }

        /************************************************
         *                  OPERATIONS                  *
         ************************************************/

        // should call this covers?
        bool contains(const Range<T> &o) const {
            return true;
        }

        bool touches(const Range<T> &o) const  {
            return true;
        }

        bool overlaps(const Range<T> &o) const  {
            return true;
        }

        Range<T> intersect(const Range<T> &o) const {
            return true;
        }

        Range<T> operator*(const Range<T> &o) const {
            return intersect(o);
        }

        bool operator<(const Range<T> o) const {
            if (lVal == o.lVal) return rVal < o.rVal;
            return lVal < o.lVal;
        }

        bool operator>(const Range<T> o) const {
            if (lVal == o.lVal) return rVal > o.rVal;
            return lVal > o.lVal;
        }

        bool operator<=(const Range<T> o) const {
            if (lVal == o.lVal) return rVal <= o.rVal;
            return lVal <= o.lVal;
        }

        bool operator>=(const Range<T> o) const {
            if (lVal == o.lVal) return rVal >= o.rVal;
            return lVal >= o.lVal;
        }

        /************************************************
         *               TRANSFORMATIONS                *
         ************************************************/
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