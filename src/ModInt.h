#ifndef MODINT_H
#define MODINT_H
#include "Constants.h"
#include "Util.h"

namespace DS {
class ModIntIndicator : MyIntegralIndicator {};

template<typename T>
using is_modint =
    typename std::conditional_t<std::is_base_of<ModIntIndicator, T>::value, std::true_type, std::false_type>;
template<typename T>
inline constexpr bool is_modint_v = is_modint<T>::value;

// @TODO Montgomery
template<uintmax_t N>
class ModInt : ModIntIndicator {
    /************************************************
     *                INITIALISATION                *
     ************************************************/

   private:
    uintmax_t val;

   public:
    // O(1)
    // if not given an initialiser value, initalise to 0
    ModInt() : val(0) {}

    // O(1)
    // Initialises a ModInt from a signed integer type
    template<typename T, std::enable_if_t<is_signed_int_v<T>, bool> = true>
    ModInt(T _num) {
        val = ((intmax_t)_num % (intmax_t)N) + N;
        val %= N;
    }

    // O(1)
    // Initialises a ModInt from an unsigned integer type
    template<typename T, std::enable_if_t<is_unsigned_int_v<T>, bool> = true>
    ModInt(T _num) : val(_num % N) {}

    /************************************************
     *                    DISPLAY                   *
     ************************************************/

    // O(1)
    // @param `out` The string representation of the graph is piped to this output stream
    // @param `newLine` Indicates whether to end with a trailing `\\n`
    friend std::ostream &operator<<(std::ostream &out, const ModInt<N> &modInt) {
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
    // @returns N
    static constexpr uintmax_t mod() {
        return N;
    }

    /************************************************
     *                   OPERATIONS                 *
     ************************************************/

    // O(1)
    // Equality
    // note: implicit conversions to modint happening here
    // e.g. ModInt<N>(x) == x + N
    friend bool operator==(const ModInt<N> &a, const ModInt<N> &b) {
        return a.num() == b.num();
    }

    // O(1)
    // Inequality
    friend bool operator!=(const ModInt<N> &a, const ModInt<N> &b) {
        return a.num() != b.num();
    }

    // O(1)
    // Standard modular addition
    friend ModInt<N> operator+(const ModInt<N> &a, const ModInt<N> &b) {
        return ModInt<N>(a.val + b.val);
    }

    // O(1)
    // Addition assignment
    friend ModInt<N> operator+=(ModInt<N> &a, const ModInt<N> b) {
        a.val += b.val;
        if (a.val >= N) a.val -= N;
        return a;
    }

    // O(1)
    // Do nothing / make copy
    ModInt<N> operator+() const {
        return ModInt<N>(val);
    }

    // O(1)
    // Prefix increment
    ModInt<N> operator++() {
        *this = *this + 1;
        return *this;
    }

    // O(1)
    // Postfix increment
    ModInt<N> operator++(int) {
        ModInt<N> tmp = *this;
        *this = *this + 1;
        return tmp;
    }

    // O(1)
    // Standard modular subtraction
    friend ModInt<N> operator-(const ModInt<N> &a, const ModInt<N> &b) {
        // some hackery to make sure our unsigned values dont overflow/underflow
        if (a.val >= b.val) {
            return ModInt<N>(a.val - b.val);
        } else {
            return -ModInt<N>(b.val - a.val);
        }
    }

    // O(1)
    // Subtraction assignment
    friend ModInt<N> operator-=(ModInt<N> &a, const ModInt<N> &b) {
        if (a.val < b.val) {
            a.val += N - b.val;
        } else {
            a.val -= b.val;
        }
        return a;
    }

    // O(1)
    // Negation
    friend ModInt<N> operator-(const ModInt<N> &a) {
        return ModInt<N>(N - a.val);
    }

    // O(1) Prefix decrement
    ModInt<N> operator--() {
        return *this = *this - 1;
    }

    // O(1) Postfix decrement
    ModInt<N> operator--(int) {
        ModInt<N> tmp = *this;
        *this = *this - 1;
        return tmp;
    }

    // O(1) Standard modular multiplication
    friend ModInt<N> operator*(const ModInt<N> a, const ModInt<N> b) {
        return ModInt<N>(a.val * b.val);
    }

    // O(1)
    // Multiplcation assignment
    friend ModInt<N> operator*=(ModInt<N> &a, const ModInt<N> &b) {
        return a = a * b;
    }

    // O(log N) Division using modular inverse
    friend ModInt<N> operator/(const ModInt<N> &a, const ModInt<N> &b) {
        return a * b.inv();
    }

    // O(log N)
    // Division assignment
    friend ModInt<N> operator/=(ModInt<N> &a, const ModInt<N> &b) {
        return a = a / b;
    }

    // O(1) Standard modulo
    // template<uintmax_t newMod> ModInt<newMod> operator%(uintmax_t o) const {
    //     return ModInt<gcd(N, o)>(val);
    // }

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

   public:
    // O(log N)
    // @returns Inverse of this number modulo `N`
    ModInt<N> inv() const {
        intmax_t x, y;
        uintmax_t g = extendedEuclidean(val, N, x, y);
        assert(g == 1 && "No solution!");
        return ModInt<N>(x);
    }

    // O(log N + gcd(denom, N))
    // solve linear congruence
    // @returns the gcd(denom, N) modints, which are solutions to val*(denom)^-1
    std::vector<ModInt<N>> allCongreunt(const ModInt<N> &denom) const {
        intmax_t u, v;
        uintmax_t d = extendedEuclidean(denom.val, N, u, v);

        std::vector<ModInt<N>> ans;
        if (val % d == 0ll) {
            ModInt<N> x0 = u * (intmax_t)(val / d);
            for (uintmax_t i = 0; i < d; i++) {
                ans.push_back(x0 + i * (N / d));
            }
        }
        return ans;
    }
};

using ModInt998244353 = ModInt<998244353>;
using ModInt1000000007 = ModInt<1000000007>;
};  // namespace DS

#endif
