#ifndef BIT_H
#define BIT_H
#include <assert.h>
#include <iostream>
#include <vector>
#include "Constants.h"
#include "Util.h"

namespace DS {
    // T() is left and right identity
    // T operator+(T, T) is associative
    // T operator-(T, T) is inverse
    template<typename T, typename Derived, size_t... Ns>
    class BITInterface_ {
private:
        T val = T();

public:
        /************************************************
         *                    DISPLAY                   *
         ************************************************/
        
        void print_helper(std::ostream& out, std::vector<T> values, size_t pos) const {
            out << values[pos];
        }

        friend std::ostream& operator<<(std::ostream& out, const BITInterface_<T, Derived, Ns...> &b) {
            out << b.val;
            return out;
        }

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/
        
        std::vector<T> data() const {
            return {val};
        }

        size_t size() const {
            return 1;
        }

        size_t dimensions() const {
            return 0;
        }

        std::vector<size_t> shape() const {
            return {};
        }

        std::vector<std::vector<size_t>> index() const {
            return {{}};
        }

        /************************************************
         *               QUREIES & UPDATES              *
         ************************************************/

        void setIndex(T v) {
            val = v;
        }
        
        void addIndex(T v) {
            val += v;
        }

        T queryIndex() const {
            return val;
        }

        void addRange(T v) {
            val += v;
        }

        T querySum() const {
            return val;
        }
    };

    template<typename T, typename Derived, size_t N, size_t... Ns>
    class BITInterface_<T, Derived, N, Ns...> {
private:
        BITInterface_<T, Derived, Ns...> bit[N + 1];

public:
        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        void print_helper(std::ostream& out, std::vector<T> values, size_t pos) const {
            out << "[ ";
            for (size_t i = 0; i < N; ++i) {
                bit[i].print_helper(out, values, pos + i * bit[i].size());
                out << ' ';
            }
            out << ']';
        }

        friend std::ostream& operator<<(std::ostream& out, const BITInterface_<T, Derived, N, Ns...> &b) {
            b.print_helper(out, static_cast<const Derived*>(&b)->data(), 0);
            return out;
        }

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        // @note derived classes should implement this
        std::vector<T> data() const;

        // O(D)
        size_t size() const {
            return N * bit[0].size();
        }

        // O(D)
        size_t dimensions() const {
            return 1 + bit[0].dimensions();
        }

        // O(D)
        std::vector<size_t> shape() const {
            std::vector<size_t> inner = bit[0].shape();
            inner.insert(inner.begin(), N);
            return inner;
        }

        // O(N^D)
        std::vector<std::vector<size_t>> index() const {
            std::vector<std::vector<size_t>> output;
            for (size_t i = 0; i < N; ++i) {
                std::vector<std::vector<size_t>> sub = bit[i].index();
                for (std::vector<size_t> subindex : sub) {
                    subindex.insert(subindex.begin(), i);
                    output.push_back(subindex);
                }
            }
            return output;
        }

        /************************************************
         *               QUREIES & UPDATES              *
         ************************************************/

        // O((2 log N)^D)
        template<typename... Args>
        void setIndex(T v, size_t pos, Args... args) {
            T currVal = static_cast<Derived*>(this)->queryIndex(pos, args...);
            static_cast<Derived*>(this)->addIndex(v - currVal, pos, args...);
        }

        // @note derived classes should implement this
        template<typename... Args>
        void addIndex(T v, size_t pos, Args... args);

        // @note derived classes should implement this
        template<typename... Args>
        T queryIndex(size_t pos, Args... args) const;
    };

    // T() is left and right identity
    // T operator+(T, T) is associative
    // T operator-(T, T) is inverse
    template<typename T, size_t... Ns>
    class BIT_RQPU: public BITInterface_<T, BIT_RQPU<T, Ns...>, Ns...> {};

    template<typename T, size_t N, size_t... Ns>
    class BIT_RQPU<T, N, Ns...>: public BITInterface_<T, BIT_RQPU<T, N, Ns...>, N, Ns...> {
private:
        BIT_RQPU<T, Ns...> bit[N + 1];

public:

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        // O((2 N log N)^D)
        std::vector<T> data() const {
            std::vector<T> output;
            for (size_t pos = 0; pos < N; ++pos) {
                std::vector<T> res(bit[0].size(), T());
                for (size_t x = pos + 1; x; x -= x & -x) res += bit[x].data();
                for (size_t x = pos + 0; x; x -= x & -x) res -= bit[x].data();
                output.insert(output.end(), res.begin(), res.end());
            }
            return output;
        }

        /************************************************
         *               QUREIES & UPDATES              *
         ************************************************/

        // O((log N)^D)
        template<typename... Args> void addIndex(T v, size_t pos, Args... args) {
            assert(0 <= pos && pos < N && "index out of range");
            for (++pos; pos <= N; pos += pos & -pos) bit[pos].addIndex(v, args...);
        }

        // O((2 log N)^D)
        template<typename... Args> T queryIndex(size_t pos, Args... args) const {
            assert(0 <= pos && pos < N && "index out of range");
            T res = T();
            for (size_t x = pos + 1; x; x -= x & -x) res += bit[x].queryIndex(args...);
            for (size_t x = pos + 0; x; x -= x & -x) res -= bit[x].queryIndex(args...);
            return res;
        }
        
        // O((2 log N)^D)
        template<typename... Args> T querySum(size_t l, size_t r, Args... args) const {
            assert(0 <= l && l <= r && r < N && "index out of range");
            T res = T();
            for (size_t x = r + 1; x; x -= x & -x) res += bit[x].querySum(args...);
            for (size_t x = l + 0; x; x -= x & -x) res -= bit[x].querySum(args...);
            return res;
        }
    };

    // T() is left and right identity
    // T operator+(T, T) is associative
    // T operator-(T, T) is inverse
    template<typename T, size_t... Ns>
    class BIT_PQRU: public BITInterface_<T, BIT_PQRU<T, Ns...>, Ns...> {};

    template<typename T, size_t N, size_t... Ns>
    class BIT_PQRU<T, N, Ns...>: public BITInterface_<T, BIT_PQRU<T, N, Ns...>, N, Ns...> {
private:
        BIT_PQRU<T, Ns...> bit[N + 1];

public:

        /************************************************
         *                  PROPERTIES                  *
         ************************************************/

        // O((N log N)^D)
        std::vector<T> data() const {
            std::vector<T> output;
            for (size_t pos = 0; pos < N; ++pos) {
                std::vector<T> res(bit[0].size(), T());
                for (size_t x = pos + 1; x <= N; x += x & -x) {
                    res += bit[x].data();
                }
                output.insert(output.end(), res.begin(), res.end());
            }
            return output;
        }

        /************************************************
         *               QUREIES & UPDATES              *
         ************************************************/

        // O((log N)^D)
        template<typename... Args> void addIndex(T v, size_t pos, Args... args) {
            assert(0 <= pos && pos < N && "index out of range");
            for (size_t x = pos + 1; x; x -= x & -x) bit[x].addIndex(v, args...);
            for (size_t x = pos + 0; x; x -= x & -x) bit[x].addIndex(-v, args...);
        }

        // O((log N)^D)
        template<typename... Args> T queryIndex(size_t pos, Args... args) const {
            assert(0 <= pos && pos < N && "index out of range");
            T res = T();
            for (size_t x = pos + 1; x <= N; x += x & -x) res += bit[x].queryIndex(args...);
            return res;
        }

        // O((2 log N)^D)
        template<typename... Args> void addRange(T v, size_t l, size_t r, Args... args) {
            assert(0 <= l && l <= r && r < N && "index out of range");
            for (size_t x = r + 1; x; x -= x & -x) bit[x].addRange(v, args...);
            for (size_t x = l + 0; x; x -= x & -x) bit[x].addRange(-v, args...);
        }
    };
};

#endif
