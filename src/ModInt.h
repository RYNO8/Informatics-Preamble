#pragma once
#include "Constants.h"
namespace DS {
    // TODO: Montgomery
    template<uint mod = MOD> class ModInt {
        /************************************************
         *                INITIALISATION                *
         ************************************************/

    private:
        ull num;

        void llInit(ll _num) {
            _num = (_num % mod) + mod;
            assert(_num >= 0);
            num = _num;
            num %= mod;
            assert(0 <= num && num < mod);
        }

    public:
        // O(1)
        // Initialises a ModInt
        ModInt(int _num = 0) {
            llInit(_num);
        }
        ModInt(uint _num = 0) : num(_num) {
            num %= mod;
            assert(0 <= num && num < mod);
        }
        ModInt(ll _num = 0) {
            llInit(_num);
        }
        ModInt(ull _num = 0) : num(_num) {
            num %= mod;
            assert(0 <= num && num < mod);
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
        const ull getNum() const {
            return num;
        }

        operator auto() const {
            return (ull)num;
        }

        /************************************************
         *                   OPERATIONS                 *
         ************************************************/
        // O(1)
        // Checks whether 2 numbers are identical
        bool operator==(ModInt o) const {
            return num == o.num;
        }

        // O(1)
        // Checks whether 2 numbers are identical
        bool operator==(ll o) const {
            return num == o;
        }

        // O(1)
        // Checks whether 2 numbers are not identical
        bool operator!=(ModInt o) const {
            return num != o.num;
        }

        // O(1)
        // Checks whether 2 numbers are not identical
        bool operator!=(ll o) const {
            return num != o;
        }

        // O(1)
        // Standard modular addition
        ModInt operator+(ModInt o) const {
            return ModInt(num + o.num);
        }

        // O(1)
        // Addition assignment
        ModInt operator+=(ModInt o) {
            *this = *this + o;
            return *this;
        }

        // O(1) Prefix increment
        ModInt operator++() {
            *this = *this + 1;
            return *this;
        }

        // O(1) Postfix increment
        ModInt operator++(int) {
            ModInt tmp = *this;
            *this = *this + 1;
            return tmp;
        }

        // O(1)
        // Standard modular subtraction
        ModInt operator-(ModInt o) const {
            return ModInt(num - o.num);
        }

        // O(1)
        // Subtraction assignment
        ModInt operator-=(ModInt o) {
            *this = *this - o;
            return *this;
        }

        // O(1) Prefix decrement
        ModInt operator--() {
            *this = *this - 1;
            return *this;
        }

        // O(1) Postfix decrement
        ModInt operator--(int) {
            ModInt tmp = *this;
            *this = *this - 1;
            return tmp;
        }

        // O(1) Standard modular multiplication
        ModInt operator*(ModInt o) const {
            return ModInt(num * o.num);
        }

        // O(1)
        // Multiplcation assignment
        void operator*=(ModInt o) {
            *this = *this * o;
        }

        // O(log x)
        // @returns Fast exponentiation of `this` to the power of `x`
        ModInt pow(ll x) const {
            assert(x >= 0);
            ModInt output(1);
            for (ModInt curr(num); x; x /= 2) {
                if (x & 1) output *= curr;
                curr *= curr;
            }
            return output;
        }

        // O(log mod)
        // @returns Inverse of this number modulo `mod`
        // ModInt modInv() const {
        //     ModInt x0(0), x1(1);
        //     ll r0 = mod, r1 = num;
        //     while (r1) {
        //         ll q = r0 / r1;
        //         x0 -= x1 * q; std::swap(x0, x1);
        //         r0 -= r1 * q; std::swap(r0, r1);
        //     }
        //     return x0;
        // }
        // TODO: make nice
        ModInt modInv() {
            if (num == 0ll) assert(false && "no inverse");
            ll A = num;
            ll M = mod;
            ll m0 = M;
            ll y = 0, x = 1;
        
            if (M == 1)
                return 0;
        
            while (A > 1) {
                if (M == 0) {
                    assert(false && "no inverse");
                }

                // q is quotient
                ll q = A / M;
                ll t = M;
        
                // m is remainder now, process same as
                // Euclid's algo
                M = A % M, A = t;
                t = y;
                
                // Update y and x
                y = x - q * y;
                x = t;
            }
        
            // Make x positive
            if (x < 0)
                x += m0;
        
            return ModInt(x);
        }

        // TODO: bezouts identity solver (i.e. congruence solver)
    };

    /************************************************
     *                    DISPLAY                   *
     ************************************************/
    template<uint mod> std::ostream& operator<<(std::ostream& out, const ModInt<mod> val) {
        val.print(out);
        return out;
    }
};
