#pragma once

// ============================================================================
// Includes
// ============================================================================

#include <vector>
#include <string>
#include <stdexcept>

// ============================================================================
// Matrix Class Declaration
// ============================================================================

template <typename T>
class Matrix {
    public:
        // Constructors and Destructors
        Matrix(size_t rows, size_t cols, T value = T()): rows_(rows), cols_(cols), data_(rows, std::vector<T>(cols, value)) {}
        ~Matrix() {}
        
        // Accessors
        size_t getRows() const {return rows_;}
        size_t getCols() const {return cols_;}

        // Element Access
        T& operator()(size_t row, size_t col);
        const T& operator()(size_t row, size_t col) const;

        // Matrix Operations
        Matrix<T> operator+(const Matrix<T>& other) const {return add(other);}
        Matrix<T> operator-(const Matrix<T>& other) const {return subtract(other);}
        Matrix<T> operator*(const Matrix<T>& other) const {return multiply(other);}
        Matrix<T> operator/(const Matrix<T>& other) const {return divide(other);}
        Matrix<T> transpose() const;
        Matrix<T> inverse() const;
        T determinant() const;
        T trace() const;
        size_t rank() const;
        size_t nullity() const;
        Matrix<T> kernel() const; // todo
        Matrix<T> image() const;// todo 
        Matrix<T> eigenvalues() const; // todo 
        Matrix<T> eigenvectors() const; // todo 
        Matrix<T> singularValues() const; // todo 

        // Matrix Properties
        bool isSquare() const {return rows_ == cols_;}
        bool isSymmetric() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isRowEchelonForm() const {return !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != 0;});}
        bool isReducedRowEchelonForm() const {return this->isRowEchelonForm() && this->isSquare();}
        bool isInvertible() const {return this->isSquare() && this->isNonSingular();}
        bool isSingular() const {return this->isSquare() && !this->isInvertible();}
        bool isNonSingular() const {return this->isSquare() && !this->isInvertible();}
        bool isPositiveDefinite() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isNegativeDefinite() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isPositiveSemidefinite() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isNegativeSemidefinite() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isIndefinite() const {return this->isSquare() && !this->isInvertible();}
        bool isHermitian() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isSkewHermitian() const {return this->isSquare() && iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isUnitary() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isOrthogonal() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        bool isIsometric() const {return this->isSquare() && !iterCheck([this](size_t i, size_t j) -> bool {return data_[i][j] != data_[j][i];});}
        
    private:
        // Member Variables
        std::vector<std::vector<T>> data_;
        size_t rows_;
        size_t cols_;

        // Internal Operations
        Matrix<T> add(const Matrix<T>& other) const;
        Matrix<T> subtract(const Matrix<T>& other) const;
        Matrix<T> multiply(const Matrix<T>& other) const;
        Matrix<T> divide(const Matrix<T>& other) const;

        // Helper Functions
        bool hasSameDimensions(const Matrix<T>& other) const {return rows_ == other.rows_ && cols_ == other.cols_;}
        template<typename Func>
        bool iterCheck(Func func) const;
        template<typename Func>
        void errorCheck(Func func, size_t row, size_t col, std::string errorMessage) const;
        template<typename Func>
        void errorCheck(Func func, std::string errorMessage) const;
        template<typename Func>
        Matrix<T> iterApply(const Matrix<T>& other, Func func) const;
        Matrix<T> iterMultiply(const Matrix<T>& other) const;
        
        // Matrix computation helpers
        void requireSquare(const std::string& operation) const;
        Matrix<T> getMinor(size_t row, size_t col) const;
        T getCofactor(size_t row, size_t col) const;
        Matrix<T> scaleByScalar(T scalar) const;
};

// ============================================================================
// Template Implementations
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Template Functions
// ----------------------------------------------------------------------------
template <typename T>
template<typename Func>
bool Matrix<T>::iterCheck(Func func) const {
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            if(func(i,j)) return true;
        }
    }
    return false;
}

template <typename T>
template<typename Func>
void Matrix<T>::errorCheck(Func func, size_t row, size_t col, std::string errorMessage) const {
    if(func(row, col)) {
        throw std::runtime_error(errorMessage);
    }
}

template <typename T>
template<typename Func>
void Matrix<T>::errorCheck(Func func, std::string errorMessage) const {
    if(func()) {
        throw std::invalid_argument(errorMessage);
    }
}

template <typename T>
template<typename Func>
Matrix<T> Matrix<T>::iterApply(const Matrix<T>& other, Func func) const {
    Matrix<T> result(rows_, cols_);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            result(i, j) = func(data_[i][j], other.data_[i][j]);
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::iterMultiply(const Matrix<T>& other) const {
    Matrix<T> result(rows_, other.cols_);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < other.cols_; j++) {
            for (size_t k = 0; k < cols_; k++) {
                result(i, j) += data_[i][k] * other.data_[k][j];
            }
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
// Element Access
// ----------------------------------------------------------------------------

template <typename T>
T& Matrix<T>::operator()(size_t row, size_t col) {
    errorCheck([this](size_t row, size_t col) -> bool {return row >= rows_ || col >= cols_;}, row, col, "Index out of range");
    return data_[row][col];
}

template <typename T>
const T& Matrix<T>::operator()(size_t row, size_t col) const {
    errorCheck([this](size_t row, size_t col) -> bool {return row >= rows_ || col >= cols_;}, row, col, "Index out of range");
    return data_[row][col];
}

// ----------------------------------------------------------------------------
// Internal Operations
// ----------------------------------------------------------------------------

template <typename T>
Matrix<T> Matrix<T>::add(const Matrix<T>& other) const {
    errorCheck([this, &other]() -> bool {return !this->hasSameDimensions(other);}, "Matrices must have the same dimensions");
    return iterApply(other, [](const T& a, const T& b) -> T {return a + b;});
}

template <typename T>
Matrix<T> Matrix<T>::subtract(const Matrix<T>& other) const {
    errorCheck([this, &other]() -> bool {return !this->hasSameDimensions(other);}, "Matrices must have the same dimensions");
    return iterApply(other, [](const T& a, const T& b) -> T {return a - b;});
}

template <typename T>
Matrix<T> Matrix<T>::multiply(const Matrix<T>& other) const {
    errorCheck([this, &other]() -> bool {return cols_ != other.rows_;}, "Matrices must have compatible dimensions for multiplication");
    return iterMultiply(other);
}

template <typename T>
Matrix<T> Matrix<T>::divide(const Matrix<T>& other) const {
    errorCheck([this, &other]() -> bool {return !this->hasSameDimensions(other);}, "Matrices must have the same dimensions");
    return iterApply(other, [](const T& a, const T& b) -> T {return a / b;});
}

// ============================================================================
// Matrix Operations
// ============================================================================

template <typename T>
Matrix<T> Matrix<T>::transpose() const {
    Matrix<T> result(cols_, rows_);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            result(j, i) = data_[i][j];
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
// Matrix Computation Helpers
// ----------------------------------------------------------------------------

template <typename T>
void Matrix<T>::requireSquare(const std::string& operation) const {
    errorCheck([this]() -> bool {return !this->isSquare();}, "Matrix must be square to " + operation);
}

template <typename T>
Matrix<T> Matrix<T>::getMinor(size_t row, size_t col) const {
    Matrix<T> minor(rows_ - 1, cols_ - 1);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            if (i != row && j != col) {
                size_t minor_i = (i < row) ? i : i - 1;
                size_t minor_j = (j < col) ? j : j - 1;
                minor(minor_i, minor_j) = data_[i][j];
            }
        }
    }
    return minor;
}

template <typename T>
T Matrix<T>::getCofactor(size_t row, size_t col) const {
    Matrix<T> minor = getMinor(row, col);
    T sign = ((row + col) % 2 == 0) ? T(1) : T(-1);
    return sign * minor.determinant();
}

template <typename T>
Matrix<T> Matrix<T>::scaleByScalar(T scalar) const {
    Matrix<T> result = *this;
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            result(i, j) *= scalar;
        }
    }
    return result;
}

// ----------------------------------------------------------------------------
// Matrix Operations
// ----------------------------------------------------------------------------

template <typename T>
T Matrix<T>::trace() const {
    requireSquare("compute trace");
    T sum = T();
    for (size_t i = 0; i < rows_; i++) {
        sum += data_[i][i];
    }
    return sum;
}

template <typename T>
T Matrix<T>::determinant() const {
    requireSquare("compute determinant");
    
    if (rows_ == 1) {
        return data_[0][0];
    }
    if (rows_ == 2) {
        return data_[0][0] * data_[1][1] - data_[0][1] * data_[1][0];
    }
    
    // For larger matrices, use recursive expansion along first row
    T det = T();
    for (size_t j = 0; j < cols_; j++) {
        T cofactor = getCofactor(0, j);
        det += data_[0][j] * cofactor;
    }
    return det;
}

template <typename T>
Matrix<T> Matrix<T>::inverse() const {
    requireSquare("compute inverse");
    
    T det = this->determinant();
    errorCheck([&det]() -> bool {return det == T(0);}, "Matrix is singular and cannot be inverted");
    
    if (rows_ == 1) {
        Matrix<T> result(1, 1);
        result(0, 0) = T(1) / data_[0][0];
        return result;
    }
    
    if (rows_ == 2) {
        Matrix<T> result(2, 2);
        T invDet = T(1) / det;
        result(0, 0) = data_[1][1] * invDet;
        result(0, 1) = -data_[0][1] * invDet;
        result(1, 0) = -data_[1][0] * invDet;
        result(1, 1) = data_[0][0] * invDet;
        return result;
    }
    
    // For larger matrices, use adjugate method
    Matrix<T> adjugate(rows_, cols_);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            adjugate(j, i) = getCofactor(i, j); // Transpose for adjugate
        }
    }
    // Multiply adjugate by 1/det (scalar multiplication)
    return adjugate.scaleByScalar(T(1) / det);
}

template <typename T>
size_t Matrix<T>::rank() const {
    // Simple rank calculation using row reduction
    // Create a copy to avoid modifying the original
    Matrix<T> temp = *this;
    size_t rank = 0;
    size_t minDim = (rows_ < cols_) ? rows_ : cols_;
    
    for (size_t col = 0; col < cols_ && rank < rows_; col++) {
        // Find pivot
        size_t pivotRow = rank;
        while (pivotRow < rows_ && temp(pivotRow, col) == T(0)) {
            pivotRow++;
        }
        
        if (pivotRow < rows_) {
            // Swap rows if needed
            if (pivotRow != rank) {
                for (size_t j = 0; j < cols_; j++) {
                    std::swap(temp(rank, j), temp(pivotRow, j));
                }
            }
            
            // Eliminate column
            for (size_t i = rank + 1; i < rows_; i++) {
                if (temp(i, col) != T(0)) {
                    T factor = temp(i, col) / temp(rank, col);
                    for (size_t j = col; j < cols_; j++) {
                        temp(i, j) -= factor * temp(rank, j);
                    }
                }
            }
            rank++;
        }
    }
    return rank;
}

template <typename T>
size_t Matrix<T>::nullity() const {
    size_t r = this->rank();
    return cols_ - r;
}

template <typename T>
Matrix<T> Matrix<T>::kernel() const {
    // Kernel (null space) - returns a matrix whose columns are basis vectors
    // This is a simplified implementation
    size_t nullityDim = this->nullity();
    
    if (nullityDim == 0) {
        // Return empty matrix or zero matrix
        return Matrix<T>(cols_, 0);
    }
    
    // For now, return identity matrix scaled by nullity
    // A proper implementation would solve Ax = 0 using row reduction
    Matrix<T> result(cols_, nullityDim, T(0));
    // TODO: Implement proper null space calculation
    return result;
}
