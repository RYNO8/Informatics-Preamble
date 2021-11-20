// @note `a` and `b` should be non-negative
ll gcd(ll a, ll b) {
	if (b == 0ll) return a;
	else return gcd(b, a % b);
}

ll gcdSlow(ll x, ll y) {
	while (y != 0) {
		ll temp = y;
		y = x % temp;
		x = temp;
	}
	return x;
}

ll lcm(ll a, ll b) {
	return a * b / gcd(a, b);
}

// O(log log A)
// @returns the number of set bits in the binary represetnation of `A`
ll bitcount(ll A) {
	A = ((A & 1431655765) + ((A & 2863311530) >> (1 << 0)));
	A = ((A & 858993459) + ((A & 3435973836) >> (1 << 1)));
	A = ((A & 252645135) + ((A & 4042322160) >> (1 << 2)));
	A = ((A & 16711935) + ((A & 4278255360) >> (1 << 3)));
	A = ((A & 65535) + ((A & 4294901760) >> (1 << 4)));
	return A;
}

// O(log A)
ll bitcountSlow(ll A) {
	ll total = 0;
	for (; A; A /= 2) total += A % 2;
	return total;
}

class ModInt {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

private:
	ll num, M;

public:
	// O(1) Initialises a ModInt
	ModInt(ll _num = 0, ll _M = MOD) {
		assert(_M > 0 && "Mod must be positive");
		M = _M;
		num = (_num + M) % _M;
		assert(0 <= num && num < M);
	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	 // O(1)
	 // @param `out` The string representation of the graph is piped to this output stream
	 // @param `newLine` Indicates whether to end with a trailing `\\n`
	void print(ostream& out = cout, bool newLine = false) {
		out << num << " ( mod " << M << " )";
		if (newLine) out << '\n';
	}

	/************************************************
	 *                   PROPERTIES                 *
	 ************************************************/

	 // O(1) Gets num
	ll getNum() {
		return num;
	}

	// O(1) Gets mod
	ll getMod() {
		return M;
	}

	/************************************************
	 *                   OPERATIONS                 *
	 ************************************************/

	 // O(1) Standard modular addition
	ModInt operator+(ModInt o) const {
		assert(M == o.M && "Can't do operations on numbers with different mods");
		return ModInt(num + o.num, M);
	}

	// O(1) Standard modular subtraction
	ModInt operator-(ModInt o) const {
		assert(M == o.M && "Can't do operations on numbers with different mods");
		return ModInt(num - o.num, M);
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
		assert(M == o.M && "Can't do operations on numbers with different mods");
		return ModInt(num * o.num, M);
	}
	void operator*=(ModInt o) {
		*this = *this * o;
	}

	// O(log o) Fast exponentiation of `this` to the power of `o`
	ModInt pow(ll o) const {
		assert(o >= 0);
		ModInt output(1, M);
		for (ModInt curr(num, M); o; o /= 2) {
			if (o & 1) output *= curr;
			curr *= curr;
		}
		return output;
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
ostream& operator<<(ostream& out, ModInt val) {
	val.print(out);
	return out;
}

template<typename T> class Matrix : public Grid<T> {
	/************************************************
	 *                INITIALISATION                *
	 ************************************************/

public:
	// lift `Grid` constructor into `Matrix` constructor
	// this means all means to initialise a `Grid` can also be used to initialise a `Matrix`
	using Grid<T>::Grid;
	//using Grid<T>::rotClockwise;

	// O(max(R, C)) Initialise a NxN identity matrix
	Matrix(int N) : Matrix(N, N, 0) {
		for (int i = 0; i < N; ++i) this->setVal(i, i, 1);
	}

	// O(RC) Initialise Matrix from Grid
	Matrix(Grid<T> g) : Grid<T>(g) {

	}

	/************************************************
	 *                    DISPLAY                   *
	 ************************************************/

	void print(ostream& out = cout, bool newLine = true) {
		out << "[\n";
		for (int r = 0; r < this->getR(); ++r) {
			out << ' ';
			for (int c = 0; c < this->getC(); ++c) out << this->getVal(r, c) << ' ';
			out << '\n';
		}
		cout << "]\n";

		if (newLine) cout << '\n';
	}

	/************************************************
	 *                   OPERATIONS                 *
	 ************************************************/

	 // O(RC) Standard matrix transposition - flip along primary diagonal
	Matrix transpose() {
		return Matrix(this->flipPrimaryDiag());
	}

	// O(RC) Standard matrix addition
	Matrix operator+(Matrix o) {
		assert(this->getC() == o.getC() && this->getR() == o.getR() && "Invalid matrix shapes for multiplication");

		Matrix* output = new Matrix<T>(*this);
		for (int r = 0; r < this->getR(); ++r) {
			for (int c = 0; c < this->getC(); ++c) {
				output->incrVal(r, c, o.getVal(r, c));
			}
		}
		return *output;
	}
	void operator+=(Matrix o) {
		*this = *this + o;
	}

	// O(RC) Stamdard matrix subtraction
	Matrix operator-(Matrix o) {
		assert(this->getC() == o.getC() && this->getR() == o.getR() && "Invalid matrix shapes for multiplication");

		Matrix* output = new Matrix<T>(*this);
		for (int r = 0; r < this->getR(); ++r) {
			for (int c = 0; c < this->getC(); ++c) {
				output->incrVal(r, c, -o.getVal(r, c));
			}
		}
		return *output;
	}
	void operator-=(Matrix o) {
		*this = *this - o;
	}

	// O(RC) Scalar multiplication
	Matrix operator*(ll mul) {
		Matrix* output = new Matrix<T>(*this);
		for (int r = 0; r < this->getR(); ++r) {
			for (int c = 0; c < this->getC(); ++c) {
				output->setVal(r, c, mul * this->getVal(r, c));
			}
		}
		return *output;
	}
	void operator*=(ll mul) {
		*this = *this * mul;
	}

	// O(N^3) Naive matrix multiplication
	Matrix operator*(Matrix o) {
		assert(this->getC() == o.getR() && "Invalid matrix shapes for multiplication");

		Matrix output(this->getR(), o.getC(), 0);
		for (int r = 0; r < this->getR(); ++r) {
			for (int c = 0; c < o.getC(); ++c) {
				for (int i = 0; i < this->getC(); ++i) {
					output.incrVal(r, c, this->getVal(r, i) * o.getVal(i, c));
				}
			}
		}
		return output;
	}
	void operator*=(Matrix o) {
		*this = *this * o;
	}

	// O(N^3 log P) Fast exponentiation, since p an integer
	Matrix pow(ll p) {
		assert(this->getR() == this->getC() && "Only square matricies can be raised to a power");

		Matrix output((int)this->getR());
		for (Matrix* curr = new Matrix(*this); p; p /= 2, *curr = (*curr) * (*curr)) {
			if (p & 1) output *= *curr;
		}
		return output;
	}
};

/************************************************
 *                    DISPLAY                   *
 ************************************************/
template<typename T> ostream& operator<<(ostream& out, Matrix<T> mat) {
	mat.print(out);
	return out;
}