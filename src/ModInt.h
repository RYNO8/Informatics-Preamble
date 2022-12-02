#pragma once
#include "Constants.h"
#include "Util.h"

namespace DS {
    // TODO: Montgomery
    template<uintmax_t mod> class ModInt {
        /************************************************
         *                INITIALISATION                *
         ************************************************/

    private:
        uintmax_t num;

    public:
        // O(1)
        // Initialises a ModInt
        template<typename T> ModInt(T _num) {
            if (std::is_signed<T>::value) {
                intmax_t unsigned_num = ((intmax_t)_num % (intmax_t)mod) + mod;
                assert(unsigned_num >= 0);
                num = unsigned_num;
                num %= mod;
            } else {
                num = _num % mod;
            }
            // invariant: num is bounded
            // implicitly also checks that mod is positive
            assert(0 <= num && num < mod);
            uintmax_t largest = mod - 1;
            assert((largest * largest) / largest == largest); // mod^2 does not overflow
        }

        // O(1)
        // if not given an initialiser value, initalise to 0
        ModInt() {
            num = 0;
        }

        /************************************************
         *                    DISPLAY                   *
         ************************************************/

        // O(1)
        // @param `out` The string representation of the graph is piped to this output stream
        // @param `newLine` Indicates whether to end with a trailing `\\n`
        void print(std::ostream& out = std::cout, bool newLine = false) const {
            out << num << " ( mod " << mod << " )";
            if (newLine) out << '\n';
        }

        /************************************************
         *                   PROPERTIES                 *
         ************************************************/

        // O(1) 
        // @returns num
        const uintmax_t getNum() const {
            return num;
        }

        operator auto() const {
            return (uintmax_t)num;
        }

        /************************************************
         *                   OPERATIONS                 *
         ************************************************/

        // O(1)
        // Standard modular addition
        ModInt<mod> operator+(ModInt<mod> o) const {
            return ModInt<mod>(num + o.num);
        }

        // O(1)
        // Addition assignment
        ModInt<mod> operator+=(ModInt<mod> o) {
            *this = *this + o;
            return *this;
        }

        // O(1) Prefix increment
        ModInt<mod> operator++() {
            *this = *this + 1;
            return *this;
        }

        // O(1) Postfix increment
        ModInt<mod> operator++(int) {
            ModInt<mod> tmp = *this;
            *this = *this + 1;
            return tmp;
        }

        // O(1)
        // Standard modular subtraction
        ModInt<mod> operator-(ModInt<mod> o) const {
            return ModInt<mod>(num - o.num);
        }

        // O(1)
        // Subtraction assignment
        ModInt<mod> operator-=(ModInt<mod> o) {
            *this = *this - o;
            return *this;
        }

        // O(1) Prefix decrement
        ModInt<mod> operator--() {
            *this = *this - 1;
            return *this;
        }

        // O(1) Postfix decrement
        ModInt<mod> operator--(int) {
            ModInt<mod> tmp = *this;
            *this = *this - 1;
            return tmp;
        }

        // O(1) Standard modular multiplication
        ModInt<mod> operator*(ModInt<mod> o) const {
            return ModInt<mod>(num * o.num);
        }

        // O(1)
        // Multiplcation assignment
        void operator*=(ModInt<mod> o) {
            *this = *this * o;
        }

        // O(1) Standard modulo
        template<uintmax_t newMod> ModInt<newMod> operator%(uintmax_t &o) const {
            return ModInt<gcd((uintmax_t)mod, o)>(num);
        }

        // O(log x)
        // @returns Fast exponentiation of `this` to the power of `x`
        ModInt<mod> pow(uintmax_t x) const {
            assert(x >= 0);
            ModInt<mod> output(1);
            for (ModInt<mod> curr(num); x; x /= 2) {
                if (x & 1) output *= curr;
                curr *= curr;
            }
            return output;
        }

private:
        // O(log mod)
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
        // O(log mod)
        // @returns Inverse of this number modulo `mod`
        ModInt<mod> modInv() const {
            intmax_t x, y;
            uintmax_t g = extendedEuclidean(num, mod, x, y);
            assert(g == 1 && "No solution!");
            return ModInt<mod>(x);
        }

        // O(log mod + gcd(denom, mod))
        // solve linear congruence
        // @returns the gcd(denom, mod) modints, which are solutions to num*(denom)^-1
        std::vector<ModInt<mod>> operator/(ModInt<mod> denom) {
            intmax_t u, v;
            uintmax_t d = extendedEuclidean(denom, mod, u, v);
            std::vector<ModInt<mod>> ans;
            if ((intmax_t)this->num % d == 0ll) {
                ModInt<mod> x0 = u * ((intmax_t)this->num / (intmax_t)d);
                for (uintmax_t i = 0; i < d; i++) {
                    ans.push_back(x0 + ModInt<mod>(i * (mod / d)));
                }
            }
            return ans;
        }
    };

    /************************************************
     *                    DISPLAY                   *
     ************************************************/
    template<uintmax_t mod> std::ostream& operator<<(std::ostream& out, const ModInt<mod> val) {
        val.print(out);
        return out;
    }
};
