template<typename T> class Polynomial {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/
private:
	vector<T> coeff;

public:
	Polynomial() {

	}

	Polynomial(vector<T> _coeff) : coeff(_coeff) {
		
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	// O(n)
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = false) {
		out << coeff;
		if (newLine) out << '\n';
	}

	// O(n), but usually O(1)
	// @param `out` The string representation of the graph is piped to this output stream
	// @param `newLine` Indicates whether to end with a trailing `\\n`
	void pprint(ostream& out = cout, bool newLine = false) {
		for (int d = deg(); d >= 0; --d) cout << coeff[d] << " x^" << d << ' ';
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/
	// O(1)
	int N() {
		return coeff.size();
	}

	// O(n), but usually O(1)
	int deg() {
		return this->norm().N() - 1;
	}

	/************************************************
	 *               BASCIC OPERATIONS              *
	 ************************************************/
	// O(n), but usually O(1)
	// @returns The polynomial without the leading 0 coefficients
	Polynomial norm() {
		Polynomial* copy = new Polynomial(*this);
		while (copy->coeff.back() == 0ll) copy->coeff.pop_back();
		return *copy;
	}

	// O(n)
	// @returns P(x) MOD x^k, in other words, do coeff[k] = 0, coeff[k+1] = 0, ...
	Polynomial mod(int k) {
		Polynomial* copy = new Polynomial(*this);
		while (copy->N() > k) copy->coeff.pop_back();
		return *copy;
	}

	// O(n)
	// @returns x^n P(1/x), in other words, reverse all coefficients
	Polynomial rev() {
		Polynomial* copy = new Polynomial(*this);
		reverse(copy->coeff.begin(), copy->coeff.end());
		return *copy;
	}

	// O(n log n)
	// @returns Q(x) = 1/P(x) mod x^t
	// in other words, the unique polynomial Q(x) such that P(x)Q(x) == 1 + x^t R(x)
	// TODO: make iterative instead
	Polynomial inv(int t) {
		assert(t != 0);
		//assert(coeff.back() != 0ll); // poly only invertible when this condition holds
		if (t == 1) return Polynomial({ coeff[0].modInv() });

		Polynomial q = this->inv((t + 1) / 2);
		Polynomial one({ 1 });
		Polynomial inv = (q - ((*this) * q - one) * q).mod(t); // idk why this expression fails when simplified
		//assert((inv * (*this)).mod(t) == one);
		return inv;
	}

	// O(n)
	// evaluate the polynomial at P(x) using Horner's method
	T eval(T x) {
		T total = T(0);
		for (int i = N() - 1; i >= 0; --i) {
			total = total * x + coeff[i];
		}
		return total;
	}

	// O(n)
	// @returns Whether 2 polynomials are identical
	bool operator==(Polynomial o) {
		Polynomial deg = this->deg();
		Polynomial oDeg = o.deg();
		if (deg != oDeg) return false;
		for (int i = 0; i < deg; ++i) {
			if (coeff[i] != o.coeff[i]) return false;
		}
		return true;
	}

	/************************************************
	 *              COMPLEX OPERATIONS              *
	 ************************************************/
	// O(n log n)
	// Number Theoric Transform (in place & iterative)
	// evaluate P(omega), P(omega^2), ..., P(omega^n), where omega is a primative `MOD` root of unity
	vector<T> NTT(int n, bool inv = false) {
		assert(popcount(n) == 1);

		vector<T> a = coeff;
		a.resize(n);

		// chance bit order, so that in place transform is possible
		for (int i = 1, j = 0; i < n; ++i) {
			int bit = n / 2;
			for (; j >= bit; bit /= 2) j -= bit;
			j += bit;
			if (i < j) swap(a[i], a[j]);
		}

		// ???
		for (int len = 2; len <= n; len *= 2) {
			assert(MOD > 2);
			T wlen = 2;
			if (inv) wlen = wlen.modInv();
			for (int i = len; i < MOD_POW; i *= 2) wlen *= wlen;

			for (int i = 0; i < n; i += len) {
				T w = 1;
				for (int j = 0; j < len / 2; ++j) {
					T u = a[i + j];
					T v = a[i + j + len / 2] * w;
					a[i + j] = u + v;
					a[i + j + len / 2] = u - v;
					w *= wlen;
				}
			}
		}

		if (inv) {
			T nInv = T(n).modInv();
			for (int i = 0; i < n; ++i) a[i] *= nInv;
		}
		return a;
	}

	// O(n log n)
	// Fast fourier transform
	// evaluate P(omega), P(omega^2), ..., P(omega^n), where omega is a complex `n`th root of unity
	vector<T> FFT(int n, bool inv = false) {

	}

	vector<T> FT(int n, bool inv = false) {
		/// TODO: check the type of T
		// if T in [int, int, ModInt]
		return NTT(n, inv);
		// else
		// return FFT(n, inv);
	}

	// O(n) Polynomial addition
	Polynomial operator+(Polynomial o) {
		Polynomial copy;
		for (int i = 0; i < max(N(), o.N()); ++i) {
			copy.coeff.push_back((i < N() ? coeff[i] : 0) + (i < o.N() ? o.coeff[i] : 0));
		}
		return copy;
	}

	// O(n) Polynomial subtraction
	Polynomial operator-(Polynomial o) {
		Polynomial copy;
		for (int i = 0; i < max(N(), o.N()); ++i) {
			copy.coeff.push_back((i < N() ? coeff[i] : 0) - (i < o.N() ? o.coeff[i] : 0));
		}
		return copy;
	}

	// O(n log n) Polynomial multiplication
	Polynomial operator*(Polynomial o) {
		int n = 1;
		while (n < 2ll * max(N(), o.N()))  n *= 2;

		vector<T> aPoints = this->FT(n, false);
		vector<T> bPoints = o.FT(n, false);

		Polynomial output;
		for (int i = 0; i < n; ++i) output.coeff.push_back(aPoints[i] * bPoints[i]);
		return output.FT(n, true);
	}

	// O(n) Polynomial scaling
	Polynomial operator*(T a) {
		Polynomial* copy = new Polynomial(*this);
		for (int i = 0; i < N(); ++i) copy->coeff[i] *= a;
		return *copy;
	}

	// O(n log n) Polynomial division
	// @returns the divisor when `this` is divided by `o`
	Polynomial operator/(Polynomial o) {
		assert(N() >= o.N());
		int divisorDeg = N() - o.N() + 1;
		Polynomial dRev = (this->rev() * o.rev().inv(divisorDeg)).mod(divisorDeg);
		Polynomial d = dRev.rev();
		return d;
	}

	// O(n log n) Polynomial divison remainder
	Polynomial operator%(Polynomial o) {
		assert(N() >= o.N());
		Polynomial r = *this - *this / o * o;
		return r;
	}

	// O(n) Addition assignment
	template <typename T> Polynomial operator+=(T o) {
		*this = (*this) + o;
		return *this;
	}

	// O(n) Subtraction assignment
	template <typename T> Polynomial operator-=(T o) {
		*this = (*this) - o;
		return *this;
	}

	// O(n log n) Multiplication assignment
	template <typename T> Polynomial operator*=(T o) {
		*this = (*this) * o;
		return *this;
	}

	// O(n log n) Division assigment
	template <typename T> Polynomial operator/=(T o) {
		*this = (*this) / o;
		return *this;
	}

	// O(n log n) Modulo assignment
	template <typename T> Polynomial operator%=(T o) {
		*this = (*this) % o;
		return *this;
	}


	vector<Polynomial> v;

	/*// O(n log^2 n)
	// TODO: this TLE's :(
	void buildSubproduct(vector<Polynomial> polys) {
		int treeSize = 1;
		while (treeSize < polys.size()) treeSize *= 2;
		v = vector<Polynomial>(treeSize);

		for (int i = 0; i < polys.size(); ++i) v[i + treeSize] = polys[i];
		for (int i = polys.size(); i < treeSize; ++i) v[i + treeSize] = Polynomial({ 1ll });

		for (int node = treeSize - 1; node >= 1; --node) {
			v[node] = v[node * 2] * v[node * 2 + 1];
		}
	}

	// @note assuming subproduct tree is already built in `v`
	vector<T> fastEval(Polynomial f, int node = 1) {
		if (2 * node >= v.size()) {
			return { f.coeff[0] };
		}
		vector<T> leftVals = fastEval(f % v[node * 2], node * 2);
		vector<T> rightVals = fastEval(f % v[node * 2 + 1], node * 2 + 1);
		leftVals.insert(leftVals.end(), rightVals.begin(), rightVals.end());
		return leftVals;
	}*/
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Polynomial<T>& val) {
	val.print(out);
	return out;
}
