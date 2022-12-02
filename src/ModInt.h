#pragma once
#include "Constants.h"
#include "Util.h"

namespace DS {
    // TODO: Montgomery
    template<uintmax_t N> class ModInt {
        /************************************************
         *                INITIALISATION                *
         ************************************************/

    private:
        uintmax_t val;

    public:
        // O(1)
        // Initialises a ModInt
        template<typename T> ModInt(T _num) {
            if (std::is_signed<T>::value) {
                intmax_t unsigned_num = ((intmax_t)_num % (intmax_t)N) + N;
                assert(unsigned_num >= 0);
                val = unsigned_num;
                val %= N;
            } else {
                val = _num % N;
            }
            // invariant: val is bounded
            // implicitly also checks that N is positive
            assert(0 <= val && val < N);
            uintmax_t largest = N - 1;
            assert((largest * largest) / largest == largest); // N^2 does not overflow
        }

        // O(1)
        // if not given an initialiser value, initalise to 0
        ModInt() {
            val = 0;
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        // O(1)
        // @param `out` The string representation of the graph is piped to this output stream
        // @param `newLine` Indicates whether to end with a trailing `\\n`
        friend std::ostream& operator<<(std::ostream& out, const ModInt modInt) {
            out << modInt.num() << " ( mod " << modInt.mod() << " )";
            return out;
        }

        /************************************************
         *                   PROPERTIES                 *
         ************************************************/

        // O(1) 
        // @returns val
        const uintmax_t num() const {
            return val;
        }

        // O(1) 
        // @returns val
        operator auto() const {
            return (uintmax_t)val;
        }

        // O(1)
        // @returns N
        static constexpr uintmax_t mod() {
            return N;
        }

        /************************************************
         *                   OPERATIONS                 *
         ************************************************/

        // O(1)
        // Standard modular addition
        ModInt<N> operator+(ModInt<N> o) const {
            return ModInt<N>(val + o.val);
        }

        // O(1)
        // Addition assignment
        ModInt<N> operator+=(ModInt<N> o) {
            *this = *this + o;
            return *this;
        }

        // O(1) Prefix increment
        ModInt<N> operator++() {
            *this = *this + 1;
            return *this;
        }

        // O(1) Postfix increment
        ModInt<N> operator++(int) {
            ModInt<N> tmp = *this;
            *this = *this + 1;
            return tmp;
        }

        // O(1)
        // Standard modular subtraction
        ModInt<N> operator-(ModInt<N> o) const {
            return ModInt<N>(val - o.val);
        }

        // O(1)
        // Subtraction assignment
        ModInt<N> operator-=(ModInt<N> o) {
            *this = *this - o;
            return *this;
        }

        // O(1) Prefix decrement
        ModInt<N> operator--() {
            *this = *this - 1;
            return *this;
        }

        // O(1) Postfix decrement
        ModInt<N> operator--(int) {
            ModInt<N> tmp = *this;
            *this = *this - 1;
            return tmp;
        }

        // O(1) Standard modular multiplication
        ModInt<N> operator*(ModInt<N> o) const {
            return ModInt<N>(val * o.val);
        }

        // O(1)
        // Multiplcation assignment
        void operator*=(ModInt<N> o) {
            *this = *this * o;
        }

        // O(1) Standard modulo
        template<uintmax_t newMod> ModInt<newMod> operator%(uintmax_t &o) const {
            return ModInt<gcd((uintmax_t)N, o)>(val);
        }

        // O(log x)
        // @returns Fast exponentiation of `this` to the power of `x`
        ModInt<N> pow(intmax_t x) const {
            ModInt<N> output(1);
            ModInt<N> curr = *this;
            if (x < 0) {
                curr = curr.inv();
                x = -x;
            }
            for (; x; x /= 2) {
                if (x & 1) output *= curr;
                curr *= curr;
            }
            return output;
        }

private:
        // O(log N)
        // extended euclidean algorithm
        // @returns gcd(a, b)
        // sets x and y such that a*x + b*y = gcd(a, b)
        static uintmax_t extendedEuclidean(uintmax_t a, uintmax_t b, intmax_t &x, intmax_t &y) {
            if (b == 0) {
                x = 1;
                y = 0;
                return a;
            }
            intmax_t x1, y1;
            uintmax_t d = extendedEuclidean(b, a % b, x1, y1);
            x = y1;
            y = x1 - y1 * (a / b);
            return d;
        }

public:
        // O(log N)
        // @returns Inverse of this number modulo `N`
        ModInt<N> inv() const {
            intmax_t x, y;
            uintmax_t g = extendedEuclidean(val, N, x, y);
            assert(g == 1 && "No solution!");
            return ModInt<N>(x);
        }

        // O(log N)
        // @returns any solution to val*(denom)^-1
         ModInt<N> operator/(ModInt<N> denom) const {
            intmax_t u, v;
            uintmax_t d = extendedEuclidean(denom, N, u, v);
            assert((intmax_t)this->val % d == 0ll);
            return ModInt<N>(u * ((intmax_t)this->val / (intmax_t)d));
        }

        // O(log N + gcd(denom, N))
        // solve linear congruence
        // @returns the gcd(denom, N) modints, which are solutions to val*(denom)^-1
        std::vector<ModInt<N>> allCongreunt(ModInt<N> denom) {
            uintmax_t d = gcd((uintmax_t)denom, N);
            std::vector<ModInt<N>> ans;
            if ((intmax_t)this->val % d == 0ll) {
                ModInt<N> x0 = *this / denom;
                for (uintmax_t i = 0; i < d; i++) {
                    ans.push_back(x0 + ModInt<N>(i * (N / d)));
                }
            }
            return ans;
        }
    };
};
