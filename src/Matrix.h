#ifndef MATRIX_H
#define MATRIX_H
#include "Constants.h"
#include "Grid.h"

namespace DS {
template<typename T, std::enable_if_t<is_my_integral_v<T>, bool> = true>
class Matrix : public Grid<T> {
    /************************************************
     *                INITIALISATION                *
     ************************************************/

   public:
    // lift `Grid` constructor into `Matrix` constructor
    // this means all means to initialise a `Grid` can also be used to initialise a `Matrix`
    using Grid<T>::Grid;

    // using Grid<T>::rotClockwise;

    // O(max(R, C))
    // Initialise a NxN identity matrix
    Matrix(int N) : Matrix(N, N, 0) {
        for (int i = 0; i < N; ++i) this->setVal(i, i, 1);
    }

    // O(RC)
    // Initialise Matrix from Grid
    Matrix(Grid<T> g) : Grid<T>(g) {}

    /************************************************
     *                    DISPLAY                   *
     ************************************************/

    // Displays the matrx
    // @param `out` The string representation of the graph is piped to this output stream
    friend std::ostream &operator<<(std::ostream &out, const Matrix<T> &matrix) {
        out << "[\n";
        for (int r = 0; r < matrix.getR(); ++r) {
            out << ' ';
            for (int c = 0; c < matrix.getC(); ++c) out << matrix.getVal(r, c) << ' ';
            out << '\n';
        }
        out << "]\n";
        return out;
    }

    /************************************************
     *                   OPERATIONS                 *
     ************************************************/

    // O(RC)
    // Standard matrix transposition - flip along primary diagonal
    Matrix transpose() const {
        return Matrix(this->flipPrimaryDiag());
    }

    // O(N)
    // @returns Matrix trace - sum of entries along main diagonal
    T trace() {
        assert(this->getR() == this->getC() && "Can't calculate trace of non-square matrix");
        T total = T(0);
        for (int i = 0; i < this->getR(); ++i) total += this->getVal(i, i);
        return total;
    }

    // O(RC)
    // Standard matrix addition
    Matrix operator+(Matrix o) const {
        assert(
            this->getC() == o.getC() && this->getR() == o.getR() && "Can't add matricies of different shapes"
        );

        Matrix *output = new Matrix<T>(*this);
        for (int r = 0; r < this->getR(); ++r) {
            for (int c = 0; c < this->getC(); ++c) {
                output->incrVal(r, c, o.getVal(r, c));
            }
        }
        return *output;
    }

    // O(RC)
    // Addition assignment
    void operator+=(Matrix o) {
        *this = *this + o;
    }

    // O(RC)
    // Stamdard matrix subtraction
    Matrix operator-(Matrix o) const {
        assert(
            this->getC() == o.getC() && this->getR() == o.getR() &&
            "Can't subtract matricies of different shapes"
        );

        Matrix *output = new Matrix<T>(*this);
        for (int r = 0; r < this->getR(); ++r) {
            for (int c = 0; c < this->getC(); ++c) {
                output->incrVal(r, c, -o.getVal(r, c));
            }
        }
        return *output;
    }

    // O(RC)
    // Subtraction assignment
    void operator-=(Matrix o) {
        *this = *this - o;
    }

    // O(RC)
    // Scalar multiplication
    Matrix operator*(ll mul) const {
        Matrix *output = new Matrix<T>(*this);
        for (int r = 0; r < this->getR(); ++r) {
            for (int c = 0; c < this->getC(); ++c) {
                output->setVal(r, c, mul * this->getVal(r, c));
            }
        }
        return *output;
    }

    // O(RC)
    // Multiplication assignment
    void operator*=(ll mul) {
        *this = *this * mul;
    }

    // O(N^3)
    // Naive matrix multiplication
    Matrix operator*(Matrix o) const {
        assert(this->getC() == o.getR() && "Invalid shapes for multiplication");

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

    // O(N^3)
    // Multiplication assignment
    void operator*=(Matrix o) {
        *this = *this * o;
    }

    // O(N^3 log P)
    // Fast exponentiation, since p an integer
    Matrix pow(ll p) const {
        assert(this->getR() == this->getC() && "Only square matricies can be raised to a power");
        assert(p >= 0 && "Inversion of matricies hasn't been implmeneted");

        Matrix output((int)this->getR());
        for (Matrix *curr = new Matrix(*this); p; p /= 2, *curr = (*curr) * (*curr)) {
            if (p & 1) output *= *curr;
        }
        return output;
    }

    // @TODO rref
    // @TODO det
    // @TODO matrix inverse
    // @TODO eigenvectors and eigenvalues
    // @TODO decompositions
};
};  // namespace DS

#endif
