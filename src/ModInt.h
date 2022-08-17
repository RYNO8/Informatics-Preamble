#pragma once

class ModInt {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

private:
	ll num;

public:
	// O(1) (if _num isn't too negative)
	// Initialises a ModInt
	ModInt(ll _num = 0) {
		while (_num < 0) _num += MOD;
		num = (_num + MOD) % MOD;
		assert(0 <= num && num < MOD);
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	// O(1)
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = false) const {
		out << num << " ( mod " << MOD << " )";
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/

	// O(1) 
	// @returns num
	const ll getNum() const {
		return num;
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

	// O(log MOD)
	// @returns Inverse of this number modulo `MOD`
	ModInt modInv() const {
		ModInt x0(0), x1(1);
		ll r0 = MOD, r1 = num;
		while (r1) {
			ll q = r0 / r1;
			x0 -= x1 * q; swap(x0, x1);
			r0 -= r1 * q; swap(r0, r1);
		}
		return x0;
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
ostream& operator<<(ostream& out, const ModInt val) {
	val.print(out);
	return out;
}

