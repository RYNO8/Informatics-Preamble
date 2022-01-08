class ModInt {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

private:
	ll num;

public:
	// O(1) Initialises a ModInt
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
	void print(ostream& out = cout, bool newLine = false) {
		out << num << " ( mod " << MOD << " )";
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/

	// O(1) Gets num
	ll getNum() {
		return num;
	}

	/************************************************
	 *                   OPERATIONS                 *
	 ************************************************/
	bool operator==(ModInt o) const {
		return num == o.num;
	}

	bool operator==(ll o) const {
		return num == o;
	}

	bool operator!=(ModInt o) const {
		return num != o.num;
	}

	bool operator!=(ll o) const {
		return num != o;
	}

	 // O(1) Standard modular addition
	ModInt operator+(ModInt o) const {
		return ModInt(num + o.num);
	}

	// O(1) Standard modular subtraction
	ModInt operator-(ModInt o) const {
		return ModInt(num - o.num);
	}

	// O(1) Addition assignment
	ModInt operator+=(ModInt o) {
		*this = *this + o;
		return *this;
	}

	// O(1) Subtraction assignment
	ModInt operator-=(ModInt o) {
		*this = *this - o;
		return *this;
	}

	// O(1) Prefix increment
	ModInt operator++() {
		*this = *this + 1;
		return *this;
	}

	// O(1) Prefix decrement
	ModInt operator--() {
		*this = *this - 1;
		return *this;
	}

	// O(1) Postfix increment
	ModInt operator++(int) {
		ModInt tmp = *this;
		*this = *this + 1;
		return tmp;
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
	void operator*=(ModInt o) {
		*this = *this * o;
	}

	// O(log o) Fast exponentiation of `this` to the power of `o`
	ModInt pow(ll o) const {
		assert(o >= 0);
		ModInt output(1);
		for (ModInt curr(num); o; o /= 2) {
			if (o & 1) output *= curr;
			curr *= curr;
		}
		return output;
	}

	ModInt modInv() {
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
ostream& operator<<(ostream& out, ModInt val) {
	val.print(out);
	return out;
}

