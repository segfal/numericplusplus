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
        Matrix<T> transpose() const; // todo
        Matrix<T> inverse() const; // todo 
        Matrix<T> determinant() const; // todo 
        Matrix<T> trace() const; // todo 
        Matrix<T> rank() const; // todo 
        Matrix<T> nullity() const; // todo 
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

