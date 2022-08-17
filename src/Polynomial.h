#pragma once

namespace DS {
	template<typename T> class Polynomial {
		/************************************************
		 *                INITIALISATION                *
		 ************************************************/
	private:
		std::vector<T> coeff;

	public:
		Polynomial() {

		}

		Polynomial(std::vector<T> _coeff) : coeff(_coeff) {
			
		}

		/************************************************
		 *                    DISPLAY                   *
		 ************************************************/

		// O(n)
		// Displays the coefficients of the polynomial
		// @param `out` The string representation of the graph is piped to this output stream
		// @param `newLine` Indicates whether to end with a trailing `\\n`
		void print(std::ostream& out = std::cout, bool newLine = false) const {
			out << coeff;
			if (newLine) out << '\n';
		}

		// O(n)
		// Displays the natural maths representation of the polynomial
		// @param `out` The string representation of the graph is piped to this output stream
		// @param `newLine` Indicates whether to end with a trailing `\\n`
		void pprint(std::ostream& out = std::cout, bool newLine = false) const {
			for (int d = deg(); d >= 0; --d) std::cout << coeff[d] << " x^" << d << ' ';
			if (newLine) out << '\n';
		}

		/************************************************
		 *                   PROPERTIES                 *
		 ************************************************/
		// O(1)
		int N() const {
			return coeff.size();
		}

		// O(n), but usually O(1)
		int deg() const {
			return this->norm().N() - 1;
		}

		/************************************************
		 *               BASCIC OPERATIONS              *
		 ************************************************/
		// O(n), but usually O(1)
		// @returns The polynomial without the leading 0 coefficients
		Polynomial<T> norm() const {
			Polynomial<T>* copy = new Polynomial<T>(*this);
			while (copy->coeff.back() == 0ll) copy->coeff.pop_back();
			return *copy;
		}

		// O(n)
		// @returns P(x) MOD x^k, in other words, do coeff[k] = 0, coeff[k+1] = 0, ...
		Polynomial<T> mod(int k) const {
			Polynomial<T>* copy = new Polynomial<T>(*this);
			while (copy->N() > k) copy->coeff.pop_back();
			return *copy;
		}

		// O(n)
		// @returns x^n P(1/x), in other words, reverse all coefficients
		Polynomial<T> rev() const {
			Polynomial<T>* copy = new Polynomial<T>(*this);
			reverse(copy->coeff.begin(), copy->coeff.end());
			return *copy;
		}

		// O(n log n)
		// @returns Q(x) = 1/P(x) mod x^t
		// in other words, the unique polynomial Q(x) such that P(x)Q(x) == 1 + x^t R(x)
		// TODO: make iterative instead
		Polynomial<T> inv(int t) const {
			assert(t != 0);
			//assert(coeff.back() != 0ll); // poly only invertible when this condition holds
			if (t == 1) return Polynomial<T>({ coeff[0].modInv() });

			Polynomial<T> q = this->inv((t + 1) / 2);
			Polynomial<T> one({ 1 });
			Polynomial<T> inv = (q - ((*this) * q - one) * q).mod(t); // idk why this expression fails when simplified
			//assert((inv * (*this)).mod(t) == one);
			return inv;
		}

		// O(n)
		// evaluate the polynomial at P(x) using Horner's method
		T eval(T x) const {
			T total = T(0);
			for (int i = N() - 1; i >= 0; --i) {
				total = total * x + coeff[i];
			}
			return total;
		}

		// O(n)
		// @returns Whether 2 polynomials are identical
		bool operator==(Polynomial<T> o) const {
			int deg = this->deg();
			int oDeg = o.deg();
			if (deg != oDeg) return false;
			for (int i = 0; i < deg; ++i) {
				if (coeff[i] != o.coeff[i]) return false;
			}
			return true;
		}

		// O(n)
		// @returns Whether 2 polynomials are not identical
		bool operator!=(Polynomial<T> o) const {
			int deg = this->deg();
			int oDeg = o.deg();
			if (deg != oDeg) return true;
			for (int i = 0; i < deg; ++i) {
				if (coeff[i] != o.coeff[i]) return true;
			}
			return false;
		}


		/************************************************
		 *              COMPLEX OPERATIONS              *
		 ************************************************/
		// O(n log n)
		// Number Theoric Transform (in place & iterative)
		// evaluate P(omega), P(omega^2), ..., P(omega^n), where omega is a primative `MOD` root of unity
		std::vector<T> NTT(int n, bool inv = false) const {
			assert(popcount(n) == 1);

			std::vector<T> a = coeff;
			a.resize(n);

			// chance bit order, so that in place transform is possible
			for (int i = 1, j = 0; i < n; ++i) {
				int bit = n / 2;
				for (; j >= bit; bit /= 2) j -= bit;
				j += bit;
				if (i < j) std::swap(a[i], a[j]);
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
		// TODO: implement
		std::vector<T> FFT(int n, bool inv = false) const {

		}

		std::vector<T> FT(int n, bool inv = false) const {
			/// TODO: check the type of T
			// if T in [int, int, ModInt]
			return NTT(n, inv);
			// else
			// return FFT(n, inv);
		}

		// O(n)
		// Polynomial<T> addition
		Polynomial<T> operator+(Polynomial<T> o) const {
			Polynomial<T> copy;
			for (int i = 0; i < std::max(N(), o.N()); ++i) {
				copy.coeff.push_back((i < N() ? coeff[i] : 0) + (i < o.N() ? o.coeff[i] : 0));
			}
			return copy;
		}

		// O(n) Addition assignment
		Polynomial<T> operator+=(Polynomial<T> o) {
			*this = (*this) + o;
			return *this;
		}

		// O(n)
		// Polynomial<T> subtraction
		Polynomial<T> operator-(Polynomial<T> o) const {
			Polynomial<T> copy;
			for (int i = 0; i < std::max(N(), o.N()); ++i) {
				copy.coeff.push_back((i < N() ? coeff[i] : 0) - (i < o.N() ? o.coeff[i] : 0));
			}
			return copy;
		}
		
		// O(n) Subtraction assignment
		Polynomial<T> operator-=(Polynomial<T> o) {
			*this = (*this) - o;
			return *this;
		}

		// O(n log n)
		// Polynomial<T> multiplication
		Polynomial<T> operator*(Polynomial<T> o) const {
			int n = 1;
			while (n < 2ll * std::max(N(), o.N()))  n *= 2;

			std::vector<T> aPoints = this->FT(n, false);
			std::vector<T> bPoints = o.FT(n, false);

			Polynomial<T> output;
			for (int i = 0; i < n; ++i) output.coeff.push_back(aPoints[i] * bPoints[i]);
			return output.FT(n, true);
		}
		
		// O(n log n) Multiplication assignment
		Polynomial<T> operator*=(Polynomial<T> o) {
			*this = (*this) * o;
			return *this;
		}

		// O(n) Polynomial<T> scaling
		Polynomial<T> operator*(T a) const {
			Polynomial* copy = new Polynomial(*this);
			for (int i = 0; i < N(); ++i) copy->coeff[i] *= a;
			return *copy;
		}
		
		// O(n) Multiplication assignment
		Polynomial<T> operator*=(T a) {
			*this = (*this) * a;
			return *this;
		}
		
		// O(n log n) Polynomial<T> division
		// @returns the divisor when `this` is divided by `o`
		Polynomial<T> operator/(Polynomial<T> o) const {
			assert(N() >= o.N());
			int divisorDeg = N() - o.N() + 1;
			Polynomial<T> dRev = (this->rev() * o.rev().inv(divisorDeg)).mod(divisorDeg);
			Polynomial<T> d = dRev.rev();
			return d;
		}

		// O(n log n) Division assigment
		Polynomial<T> operator/=(Polynomial<T> o) {
			*this = (*this) / o;
			return *this;
		}

		// O(n log n) Polynomial<T> divison remainder
		// @returns the remainer when `this` is divided by `o`
		Polynomial<T> operator%(Polynomial<T> o) const {
			assert(N() >= o.N());
			Polynomial<T> r = *this - *this / o * o;
			return r;
		}

		// O(n log n) Modulo assignment
		Polynomial<T> operator%=(Polynomial<T>o) {
			*this = (*this) % o;
			return *this;
		}

	private:
		// O(n log^2 n)
		// Build subproduct tree, where the root is v[1], left child is 2n, and right child is 2n+1
		// Each node contains the product of all the leaf node polynomials that it covers
		std::vector<Polynomial<T>> buildSubproduct(std::vector<T> points) const {
			int treeSize = 1;
			while (treeSize < points.size()) treeSize *= 2;
			std::vector<Polynomial<T>> v = std::vector<Polynomial<T>>(treeSize);

			for (int i = 0; i < points.size(); ++i) v[i + treeSize] = Polynomial<T>({ -points[i], 1ll });
			for (int i = points.size(); i < treeSize; ++i) v[i + treeSize] = Polynomial<T>({ 1ll });

			for (int node = treeSize - 1; node >= 1; --node) {
				v[node] = v[node * 2] * v[node * 2 + 1];
			}
		}

		// O(n log^2 n)
		// Divide and conquer, degree of f = degree of v[node] - 1
		std::vector<T> fastEval_(std::vector<Polynomial<T>> &v, int node = 1) const {
			if (2 * node >= v.size()) {
				return { coeff[0] };
			}
			std::vector<T> leftVals = (*this % v[node * 2]).fastEval(v, node * 2);
			std::vector<T> rightVals = (*this % v[node * 2 + 1]).fastEval(v, node * 2 + 1);
			leftVals.insert(leftVals.end(), rightVals.begin(), rightVals.end());
			return leftVals;
		}

	public:
		// O(n log^2 n)
		// Fast multipoint evaluation
		// @return P(points[0]), P(points[1]), ...
		std::vector<T> fastEval(std::vector<T> points) const {
			std::vector<Polynomial<T>> v = buildSubproduct(points);
			return fastEval_(v);
		}
	};

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/
	template<typename T> std::ostream& operator<<(std::ostream& out, const Polynomial<T>& val) {
		val.print(out);
		return out;
	}
};
