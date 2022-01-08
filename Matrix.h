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