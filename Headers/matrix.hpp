#pragma once

// ============================================================================
// Includes
// ============================================================================

#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>

// ============================================================================
// Matrix Class Declaration
// ============================================================================

template <typename T>
class Matrix {
    public:
        // Constructors and Destructors
        Matrix(size_t rows, size_t cols, T value = T()): rows_(rows), cols_(cols), data_(rows, std::vector<T>(cols, value)) {}
        Matrix(const std::vector<std::vector<T>>& data): rows_(data.size()), cols_(data[0].size()), data_(data) {}
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
        
        Matrix<T> rowReduce() const;
        Matrix<T> columnReduce() const {return this->transpose().rowReduce().transpose();}
        Matrix<T> reducedRowEchelonForm() const;
        Matrix<T> transpose() const;
        Matrix<T> inverse() const;
        T determinant() const;
        T trace() const;
        size_t rank() const;
        size_t nullity() const; // incomplete
        Matrix<T> kernel() const; // incomplete
        Matrix<T> image() const;// incomplete
        Matrix<T> eigenvalues() const; // incomplete
        Matrix<T> eigenvectors() const; // incomplete
        Matrix<T> singularValues() const; // incomplete

        // Matrix Properties
        bool isSquare() const {return rows_ == cols_;}
        bool isSymmetric() const;
        bool isRowEchelonForm() const;
        bool isReducedRowEchelonForm() const;
        bool isInvertible() const {return this->isNonSingular();}
        bool isSingular() const;
        bool isNonSingular() const;
        bool isPositiveDefinite() const;
        bool isNegativeDefinite() const;
        bool isPositiveSemidefinite() const;
        bool isNegativeSemidefinite() const;
        bool isIndefinite() const;
        bool isHermitian() const {return this->isSymmetric();}
        bool isSkewHermitian() const;
        bool isUnitary() const {return this->isOrthogonal();}
        bool isOrthogonal() const;
        bool isIsometric() const;
        
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
    // Kernel (null space) - returns a matrix whose columns are basis vectors for the null space
    // Solves Ax = 0 using RREF
    size_t nullityDim = this->nullity();
    
    if (nullityDim == 0) {
        // Return empty matrix (cols_ x 0)
        return Matrix<T>(cols_, 0);
    }
    
    // Get RREF to identify pivot columns and free variables
    Matrix<T> rref = this->reducedRowEchelonForm();
    
    // Find pivot columns (columns with leading 1s)
    std::vector<bool> isPivotCol(cols_, false);
    std::vector<size_t> pivotCols;
    size_t currentRow = 0;
    
    for (size_t col = 0; col < cols_ && currentRow < rows_; col++) {
        if (rref(currentRow, col) != T(0)) {
            // Found a pivot
            isPivotCol[col] = true;
            pivotCols.push_back(col);
            currentRow++;
        }
    }
    
    // Create basis vectors for free variables (non-pivot columns)
    Matrix<T> result(cols_, nullityDim, T(0));
    size_t basisIdx = 0;
    
    for (size_t col = 0; col < cols_; col++) {
        if (!isPivotCol[col]) {
            // This is a free variable - create a basis vector
            result(col, basisIdx) = T(1); // Set free variable to 1
            
            // Back-substitute to find values of pivot variables
            // For each row in RREF, solve for pivot variables
            for (int row = static_cast<int>(pivotCols.size()) - 1; row >= 0; row--) {
                size_t pivotCol = pivotCols[row];
                T sum = T(0);
                
                // Sum all non-pivot contributions
                for (size_t j = pivotCol + 1; j < cols_; j++) {
                    if (!isPivotCol[j] || j == col) {
                        sum += rref(row, j) * result(j, basisIdx);
                    }
                }
                
                // Solve: pivotVar = -sum (since pivot is 1 in RREF)
                result(pivotCol, basisIdx) = -sum;
            }
            
            basisIdx++;
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::image() const {
   
    size_t rankDim = this->rank();
    
    if (rankDim == 0) {
        return Matrix<T>(rows_, 0);
    }
    
    Matrix<T> rref = this->reducedRowEchelonForm();
    
    std::vector<size_t> pivotCols;
    size_t currentRow = 0;
    
    for (size_t col = 0; col < cols_ && currentRow < rows_; col++) {
        if (rref(currentRow, col) != T(0)) {
            pivotCols.push_back(col);
            currentRow++;
        }
    }
    
    Matrix<T> result(rows_, rankDim);
    for (size_t i = 0; i < pivotCols.size(); i++) {
        size_t pivotCol = pivotCols[i];
        for (size_t row = 0; row < rows_; row++) {
            result(row, i) = data_[row][pivotCol];
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::eigenvalues() const {
   
    requireSquare("compute eigenvalues");
    
    if (rows_ == 1) {
        Matrix<T> result(1, 1);
        result(0, 0) = data_[0][0];
        return result;
    }
    
    if (rows_ == 2) {
        T trace = this->trace();
        T det = this->determinant();
        
        T discriminant = trace * trace - T(4) * det;
        
        Matrix<T> result(2, 1);
        
        if (discriminant < T(0)) {
            T realPart = trace / T(2);
            T imagPart = std::sqrt(-discriminant) / T(2);
            result(0, 0) = realPart;
            result(1, 0) = realPart;
        } else {
            T sqrtDisc = std::sqrt(discriminant);
            result(0, 0) = (trace + sqrtDisc) / T(2);
            result(1, 0) = (trace - sqrtDisc) / T(2);
        }
        
        return result;
    }
    
    Matrix<T> result(rows_, 1, T(0));
    
    Matrix<T> v(rows_, 1, T(1));
    
    T norm = T(0);
    for (size_t i = 0; i < rows_; i++) {
        norm += v(i, 0) * v(i, 0);
    }
    norm = std::sqrt(norm);
    if (norm != T(0)) {
        for (size_t i = 0; i < rows_; i++) {
            v(i, 0) /= norm;
        }
    }
    
    const size_t maxIter = 100;
    const T tolerance = T(1e-10);
    
    for (size_t iter = 0; iter < maxIter; iter++) {
        Matrix<T> Av = (*this) * v;
        
        T numerator = T(0);
        T denominator = T(0);
        for (size_t i = 0; i < rows_; i++) {
            numerator += v(i, 0) * Av(i, 0);
            denominator += v(i, 0) * v(i, 0);
        }
        
        T eigenvalue = (denominator != T(0)) ? numerator / denominator : T(0);
        result(0, 0) = eigenvalue;
        
        norm = T(0);
        for (size_t i = 0; i < rows_; i++) {
            norm += Av(i, 0) * Av(i, 0);
        }
        norm = std::sqrt(norm);
        
        if (norm < tolerance) break;
        
        for (size_t i = 0; i < rows_; i++) {
            v(i, 0) = Av(i, 0) / norm;
        }
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::eigenvectors() const {
    requireSquare("compute eigenvectors");
    
    Matrix<T> eigenvals = this->eigenvalues();
    
    Matrix<T> result(rows_, rows_, T(0));
    
    if (rows_ == 1) {
        result(0, 0) = T(1);
        return result;
    }
    
    if (rows_ == 2) {
        for (size_t i = 0; i < 2 && i < eigenvals.getRows(); i++) {
            T lambda = eigenvals(i, 0);
            
            Matrix<T> A_minus_lambdaI = *this;
            A_minus_lambdaI(0, 0) -= lambda;
            A_minus_lambdaI(1, 1) -= lambda;
            
            Matrix<T> eigenvec = A_minus_lambdaI.kernel();
            
            if (eigenvec.getCols() > 0) {
                T norm = T(0);
                for (size_t j = 0; j < eigenvec.getRows(); j++) {
                    norm += eigenvec(j, 0) * eigenvec(j, 0);
                }
                norm = std::sqrt(norm);
                if (norm != T(0)) {
                    for (size_t j = 0; j < eigenvec.getRows(); j++) {
                        result(j, i) = eigenvec(j, 0) / norm;
                    }
                } else {
                    if (A_minus_lambdaI(0, 0) != T(0) || A_minus_lambdaI(1, 0) != T(0)) {
                        result(0, i) = -A_minus_lambdaI(1, 0);
                        result(1, i) = A_minus_lambdaI(0, 0);
                        norm = std::sqrt(result(0, i) * result(0, i) + result(1, i) * result(1, i));
                        if (norm != T(0)) {
                            result(0, i) /= norm;
                            result(1, i) /= norm;
                        }
                    } else {
                        result(0, i) = T(1);
                        result(1, i) = T(0);
                    }
                }
            } else {
                result(0, i) = T(1);
                result(1, i) = T(0);
            }
        }
        return result;
    }
    
    Matrix<T> v(rows_, 1, T(1));
    
    T norm = T(0);
    for (size_t i = 0; i < rows_; i++) {
        norm += v(i, 0) * v(i, 0);
    }
    norm = std::sqrt(norm);
    if (norm != T(0)) {
        for (size_t i = 0; i < rows_; i++) {
            v(i, 0) /= norm;
        }
    }
    
    const size_t maxIter = 100;
    for (size_t iter = 0; iter < maxIter; iter++) {
        Matrix<T> Av = (*this) * v;
        
        norm = T(0);
        for (size_t i = 0; i < rows_; i++) {
            norm += Av(i, 0) * Av(i, 0);
        }
        norm = std::sqrt(norm);
        
        if (norm < T(1e-10)) break;
        
        for (size_t i = 0; i < rows_; i++) {
            v(i, 0) = Av(i, 0) / norm;
        }
    }
    
    for (size_t i = 0; i < rows_; i++) {
        result(i, 0) = v(i, 0);
    }
    
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::singularValues() const {

    size_t minDim = (rows_ < cols_) ? rows_ : cols_;
    Matrix<T> result(minDim, 1);
    
    Matrix<T> ATA = this->transpose() * (*this);
    
    Matrix<T> eigenvals = ATA.eigenvalues();
    
    size_t count = (eigenvals.getRows() < minDim) ? eigenvals.getRows() : minDim;
    for (size_t i = 0; i < count; i++) {
        T eigenval = eigenvals(i, 0);
        if (eigenval < T(0)) eigenval = T(0);
        result(i, 0) = std::sqrt(eigenval);
    }
    
    for (size_t i = 0; i < count; i++) {
        for (size_t j = i + 1; j < count; j++) {
            if (result(i, 0) < result(j, 0)) {
                std::swap(result(i, 0), result(j, 0));
            }
        }
    }
    
    return result;
}
template <typename T>
Matrix<T> Matrix<T>::rowReduce() const {
    Matrix<T> result = *this;
    size_t rank = 0;
    
    for (size_t col = 0; col < cols_ && rank < rows_; col++) {
        size_t pivotRow = rank;
        while (pivotRow < rows_ && result(pivotRow, col) == T(0)) {
            pivotRow++;
        }
        
        if (pivotRow < rows_) {
            if (pivotRow != rank) {
                for (size_t j = 0; j < cols_; j++) {
                    std::swap(result(rank, j), result(pivotRow, j));
                }
            }
            
            for (size_t i = rank + 1; i < rows_; i++) {
                if (result(i, col) != T(0)) {
                    T factor = result(i, col) / result(rank, col);
                    for (size_t j = col; j < cols_; j++) {
                        result(i, j) -= factor * result(rank, j);
                    }
                }
            }
            rank++;
        }
    }
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::reducedRowEchelonForm() const {
    Matrix<T> result = this->rowReduce();
    size_t rank = 0;
    
    for (size_t i = 0; i < rows_; i++) {
        bool isNonZeroRow = false;
        for (size_t j = 0; j < cols_; j++) {
            if (result(i, j) != T(0)) {
                isNonZeroRow = true;
                break;
            }
        }
        if (isNonZeroRow) rank++;
    }
    
for (int i = static_cast<int>(rank) - 1; i >= 0; i--) {
        size_t pivotCol = 0;
        while (pivotCol < cols_ && result(i, pivotCol) == T(0)) {
            pivotCol++;
        }
        
        if (pivotCol < cols_ && result(i, pivotCol) != T(0)) {
            T pivot = result(i, pivotCol);
            for (size_t j = pivotCol; j < cols_; j++) {
                result(i, j) /= pivot;
            }
            
            for (int k = i - 1; k >= 0; k--) {
                if (result(k, pivotCol) != T(0)) {
                    T factor = result(k, pivotCol);
                    for (size_t j = pivotCol; j < cols_; j++) {
                        result(k, j) -= factor * result(i, j);
                    }
                }
            }
        }
    }
    return result;
}

// ============================================================================
// Matrix Property Implementations
// ============================================================================

template <typename T>
bool Matrix<T>::isSymmetric() const {
    if (!this->isSquare()) return false;
    return !iterCheck([this](size_t i, size_t j) -> bool {
        return data_[i][j] != data_[j][i];
    });
}

template <typename T>
bool Matrix<T>::isRowEchelonForm() const {
 
    
    size_t lastPivotCol = SIZE_MAX;
    
    for (size_t i = 0; i < rows_; i++) {
        size_t pivotCol = cols_;
        for (size_t j = 0; j < cols_; j++) {
            if (data_[i][j] != T(0)) {
                pivotCol = j;
                break;
            }
        }
        
        if (pivotCol == cols_) {
            for (size_t k = i + 1; k < rows_; k++) {
                for (size_t j = 0; j < cols_; j++) {
                    if (data_[k][j] != T(0)) return false;
                }
            }
            break;
        } else {
            if (pivotCol <= lastPivotCol && lastPivotCol != SIZE_MAX) {
                return false;
            }
            lastPivotCol = pivotCol;
            
            for (size_t k = i + 1; k < rows_; k++) {
                if (data_[k][pivotCol] != T(0)) return false;
            }
        }
    }
    
    return true;
}

template <typename T>
bool Matrix<T>::isReducedRowEchelonForm() const {
    
    if (!this->isRowEchelonForm()) return false;
    
    for (size_t i = 0; i < rows_; i++) {
        size_t pivotCol = cols_;
        for (size_t j = 0; j < cols_; j++) {
            if (data_[i][j] != T(0)) {
                pivotCol = j;
                break;
            }
        }
        
        if (pivotCol == cols_) continue;
        
        if (data_[i][pivotCol] != T(1)) return false;
        
        for (size_t k = 0; k < rows_; k++) {
            if (k != i && data_[k][pivotCol] != T(0)) return false;
        }
    }
    
    return true;
}

template <typename T>
bool Matrix<T>::isNonSingular() const {
    if (!this->isSquare()) return false;
    return this->determinant() != T(0);
}

template <typename T>
bool Matrix<T>::isSingular() const {
    if (!this->isSquare()) return false;
    return this->determinant() == T(0);
}

template <typename T>
bool Matrix<T>::isPositiveDefinite() const {
    if (!this->isSquare() || !this->isSymmetric()) return false;
    
    for (size_t k = 1; k <= rows_; k++) {
        Matrix<T> submatrix(k, k);
        for (size_t i = 0; i < k; i++) {
            for (size_t j = 0; j < k; j++) {
                submatrix(i, j) = data_[i][j];
            }
        }
        if (submatrix.determinant() <= T(0)) return false;
    }
    return true;
}

template <typename T>
bool Matrix<T>::isNegativeDefinite() const {
    if (!this->isSquare() || !this->isSymmetric()) return false;
    
    Matrix<T> negA = this->scaleByScalar(T(-1));
    return negA.isPositiveDefinite();
}

template <typename T>
bool Matrix<T>::isPositiveSemidefinite() const {
    if (!this->isSquare() || !this->isSymmetric()) return false;
    
    if (rows_ <= 2) {
        Matrix<T> eigenvals = this->eigenvalues();
        for (size_t i = 0; i < eigenvals.getRows(); i++) {
            if (eigenvals(i, 0) < T(0)) return false;
        }
        return true;
    }
    
   
    return false;   
}

template <typename T>
bool Matrix<T>::isNegativeSemidefinite() const {
    if (!this->isSquare() || !this->isSymmetric()) return false;
    
    Matrix<T> negA = this->scaleByScalar(T(-1));
    return negA.isPositiveSemidefinite();
}

template <typename T>
bool Matrix<T>::isIndefinite() const {
    if (!this->isSquare() || !this->isSymmetric()) return false;
    
    if (rows_ <= 2) {
        Matrix<T> eigenvals = this->eigenvalues();
        bool hasPositive = false;
        bool hasNegative = false;
        for (size_t i = 0; i < eigenvals.getRows(); i++) {
            if (eigenvals(i, 0) > T(0)) hasPositive = true;
            if (eigenvals(i, 0) < T(0)) hasNegative = true;
        }
        return hasPositive && hasNegative;
    }
    
    return false; 
}

template <typename T>
bool Matrix<T>::isSkewHermitian() const {
    if (!this->isSquare()) return false;
    
    return !iterCheck([this](size_t i, size_t j) -> bool {
        if (i == j) return data_[i][j] != T(0);
        return data_[i][j] != -data_[j][i];
    });
}

template <typename T>
bool Matrix<T>::isOrthogonal() const {
    if (!this->isSquare()) return false;
    
    Matrix<T> AAT = (*this) * this->transpose();
    
    const T tolerance = T(1e-10);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            T expected = (i == j) ? T(1) : T(0);
            if (std::abs(AAT(i, j) - expected) > tolerance) return false;
        }
    }
    return true;
}

template <typename T>
bool Matrix<T>::isIsometric() const {

    if (!this->isSquare()) return false;
    
    Matrix<T> ATA = this->transpose() * (*this);
    
    const T tolerance = T(1e-10);
    for (size_t i = 0; i < rows_; i++) {
        for (size_t j = 0; j < cols_; j++) {
            T expected = (i == j) ? T(1) : T(0);
            if (std::abs(ATA(i, j) - expected) > tolerance) return false;
        }
    }
    return true;
}
