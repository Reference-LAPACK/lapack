
#ifndef LAPACK_CPP_MATRIX_HPP
#define LAPACK_CPP_MATRIX_HPP

#include <algorithm>
#include <cassert>
#include <chrono>
#include <functional>
#include <iostream>
#include <random>

#include "lapack_cpp/base/types.hpp"
#include "lapack_cpp/base/enums.hpp"
#include "lapack_cpp/base/memory_block.hpp"
#include "lapack_cpp/base/vector.hpp"

namespace lapack_cpp {

// Forward declaration of ConstMatrix
template <typename T, Layout layout = Layout::ColMajor, typename idx_t = lapack_idx_t>
class ConstMatrix;

/**
 * A matrix class that can be used to store matrices of arbitrary size.
 *
 * The elements of the matrix are stored in column-major order
 *
 * @tparam T this is a template parameter that specifies the type of the
 *           elements of the vector.
 */
template <typename T, Layout layout = Layout::ColMajor, typename idx_t = lapack_idx_t>
class Matrix {
   public:
    typedef T value_type;
    typedef idx_t index_type;
    // Constructor for a matrix of size m x n
    template <bool aligned>
    Matrix(const idx_t m,
           const idx_t n,
           const MemoryBlock<T, idx_t, aligned>& m_block)
        : m_(m),
          n_(n),
          ldim_(calc_ld<T, idx_t, aligned>(layout == Layout::ColMajor ? m_ : n_)),
          data_(m_block.ptr())
    {
        assert(m_block.size() >=
               ldim_ * (layout == Layout::ColMajor ? n_ : m_));
    }

    // Constructor for a matrix of size m x n
    template <bool aligned>
    Matrix(const idx_t m,
           const idx_t n,
           const MemoryBlock<T, idx_t, aligned>& m_block,
           const idx_t ldim)
        : m_(m), n_(n), ldim_(ldim), data_(m_block.ptr())
    {
        assert(m_block.size() >=
               ldim_ * (layout == Layout::ColMajor ? n_ : m_));
    }

    // Constructor for a matrix of size m x n
    Matrix(const idx_t m, const idx_t n, T* ptr, const idx_t ldim)
        : m_(m), n_(n), ldim_(ldim), data_(ptr)
    {
        assert(ldim >= (layout == Layout::ColMajor ? m : n));
    }

    // Copy constructor
    Matrix(const Matrix& m)
        : m_(m.num_rows()), n_(m.num_columns()), ldim_(m.ldim()), data_(m.ptr())
    {}

    // Make a const promotion of the matrix
    ConstMatrix<T, layout, idx_t> as_const() const
    {
        return ConstMatrix<T, layout, idx_t>(m_, n_, data_, ldim_);
    }

    // Assignment operator
    void operator=(const Matrix& m)
    {
        assert(m.num_rows() == m_);
        assert(m.num_columns() == n_);
        for (idx_t i = 0; i < m_; ++i) {
            for (idx_t j = 0; j < n_; ++j) {
                (*this)(i,j) = m(i, j);
            }
        }
    }

    // Returns the number of rows of the matrix
    inline idx_t num_rows() const { return m_; }

    // Returns the number of columns of the matrix
    inline idx_t num_columns() const { return n_; }

    // Returns the leading dimension of the matrix
    inline idx_t ldim() const { return ldim_; }

    // Returns the size of the matrix (i.e. the number of elements m * n)
    inline idx_t size() const { return n_ * m_; }

    // Returns the ij-th element of the matrix
    inline T& operator()(idx_t i, idx_t j) const
    {
        assert(i < m_);
        assert(j < n_);
        if (layout == Layout::ColMajor)
            return data_[i + j * ldim_];
        else
            return data_[i * ldim_ + j];
    }

    /**
     * Return a submatrix of the matrix.
     *
     * @param i1 index of the first row of the submatrix
     * @param i2 index of the last row of the submatrix (not inclusive)
     * @param j1 index of the first column of the submatrix
     * @param j2 index of the last column of the submatrix (not inclusive)
     *
     * @example
     * Matrix<float> A(5,5);
     * auto B = A.submatrix(1, 4, 1, 4); // B is a 3x3 submatrix of A
     */
    inline Matrix submatrix(idx_t i1, idx_t i2, idx_t j1, idx_t j2) const
    {
        assert(i1 >= 0);
        assert(i2 <= m_);
        assert(j1 >= 0);
        assert(j2 <= n_);
        assert(i1 < i2);
        assert(j1 < j2);
        if (layout == Layout::ColMajor)
            return Matrix(i2 - i1, j2 - j1, &data_[i1 + j1 * ldim_], ldim_);
        else
            return Matrix(i2 - i1, j2 - j1, &data_[i1 * ldim_ + j1], ldim_);
    }

    /**
     * Return a column of the matrix as a vector
     *
     * @param i index of the column
     */
    inline Vector<T, idx_t> column(idx_t i) const
    {
        assert(i >= 0);
        assert(i < n_);
        if (layout == Layout::ColMajor)
            return Vector<T, idx_t>(m_, &data_[i * ldim_], 1);
        else
            return Vector<T, idx_t>(m_, &data_[i], ldim_);
    }

    /**
     * Return a row of the matrix as a vector
     *
     * @param i index of the row
     */
    inline Vector<T, idx_t> row(idx_t i) const
    {
        assert(i >= 0);
        assert(i < m_);
        if (layout == Layout::ColMajor)
            return Vector<T, idx_t>(n_, &data_[i], ldim_);
        else
            return Vector<T, idx_t>(n_, &data_[i * ldim_], 1);
    }

    /**
     * Return a pointer to the data of the matrix.
     */
    inline T* ptr() const { return data_; }

   private:
    // Number of rows of the matrix
    const idx_t m_;
    // Number of columns of the matrix
    const idx_t n_;
    // Leading dimension of the matrix
    const idx_t ldim_;
    // Pointer to the data
    T* data_;
};

template <typename T, Layout layout, typename idx_t >
class ConstMatrix {
   public:
    typedef T value_type;
    typedef idx_t index_type;
    // Constructor for a matrix of size m x n
    template <bool aligned>
    ConstMatrix(const idx_t m,
                const idx_t n,
                const MemoryBlock<T, idx_t, aligned>& m_block)
        : m_(m),
          n_(n),
          ldim_(calc_ld<T, idx_t, aligned>(layout == Layout::ColMajor ? m_ : n_)),
          data_(m_block.ptr())
    {
        assert(m_block.size() >=
               ldim_ * (layout == Layout::ColMajor ? n_ : m_));
    }

    // Constructor for a matrix of size m x n
    template <bool aligned>
    ConstMatrix(const idx_t m,
                const idx_t n,
                const MemoryBlock<T, idx_t, aligned>& m_block,
                const idx_t ldim)
        : m_(m), n_(n), ldim_(ldim), data_(m_block.ptr())
    {
        assert(m_block.size() >=
               ldim_ * (layout == Layout::ColMajor ? n_ : m_));
    }

    // Constructor for a matrix of size m x n
    ConstMatrix(const idx_t m, const idx_t n, const T* ptr, const idx_t ldim)
        : m_(m), n_(n), ldim_(ldim), data_(ptr)
    {
        assert(ldim >= (layout == Layout::ColMajor ? m : n));
    }

    // Copy constructor
    ConstMatrix(const ConstMatrix& m)
        : m_(m.num_rows()), n_(m.num_columns()), ldim_(m.ldim()), data_(m.ptr())
    {}

    // Const promotion constructor
    ConstMatrix(const Matrix<T, layout, idx_t>& m)
        : m_(m.num_rows()), n_(m.num_columns()), ldim_(m.ldim()), data_(m.ptr())
    {}

    // Returns the number of rows of the matrix
    inline idx_t num_rows() const { return m_; }

    // Returns the number of columns of the matrix
    inline idx_t num_columns() const { return n_; }

    // Returns the leading dimension of the matrix
    inline idx_t ldim() const { return ldim_; }

    // Returns the size of the matrix (i.e. the number of elements m * n)
    inline idx_t size() const { return n_ * m_; }

    // Returns the ij-th element of the matrix
    inline T operator()(idx_t i, idx_t j) const
    {
        assert(i < m_);
        assert(j < n_);
        if (layout == Layout::ColMajor)
            return data_[i + j * ldim_];
        else
            return data_[i * ldim_ + j];
    }

    /**
     * Return a submatrix of the matrix.
     *
     * @param i1 index of the first row of the submatrix
     * @param i2 index of the last row of the submatrix (not inclusive)
     * @param j1 index of the first column of the submatrix
     * @param j2 index of the last column of the submatrix (not inclusive)
     *
     * @example
     * Matrix<float> A(5,5);
     * auto B = A.submatrix(1, 4, 1, 4); // B is a 3x3 submatrix of A
     */
    inline ConstMatrix submatrix(idx_t i1, idx_t i2, idx_t j1, idx_t j2) const
    {
        assert(i1 >= 0);
        assert(i2 <= m_);
        assert(j1 >= 0);
        assert(j2 <= n_);
        assert(i1 < i2);
        assert(j1 < j2);
        if (layout == Layout::ColMajor)
            return ConstMatrix(i2 - i1, j2 - j1, &data_[i1 + j1 * ldim_],
                               ldim_);
        else
            return ConstMatrix(i2 - i1, j2 - j1, &data_[i1 * ldim_ + j1],
                               ldim_);
    }

    /**
     * Return a column of matrix as a vector
     *
     * @param i index of the column
     */
    inline ConstVector<T, idx_t> column(idx_t i) const
    {
        assert(i >= 0);
        assert(i < n_);
        if (layout == Layout::ColMajor)
            return ConstVector<T, idx_t>(m_, &data_[i * ldim_], 1);
        else
            return ConstVector<T, idx_t>(m_, &data_[i], ldim_);
    }

    /**
     * Return a row of matrix as a vector
     *
     * @param i index of the row
     */
    inline ConstVector<T, idx_t> row(idx_t i) const
    {
        assert(i >= 0);
        assert(i < m_);
        if (layout == Layout::ColMajor)
            return ConstVector<T, idx_t>(n_, &data_[i], ldim_);
        else
            return ConstVector<T, idx_t>(n_, &data_[i * ldim_], 1);
    }

    /**
     * Return a pointer to the data of the matrix.
     */
    inline const T* ptr() const { return data_; }

   private:
    // Number of rows of the matrix
    const idx_t m_;
    // Number of columns of the matrix
    const idx_t n_;
    // Leading dimension of the matrix
    const idx_t ldim_;
    // Pointer to the data
    const T* data_;
};

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_MATRIX_HPP