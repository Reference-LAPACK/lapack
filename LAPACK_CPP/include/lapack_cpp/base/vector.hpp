
#ifndef LAPACK_CPP_VECTOR_HPP
#define LAPACK_CPP_VECTOR_HPP

#include <algorithm>
#include <cassert>
#include <chrono>

#include "lapack_cpp/base/types.hpp"
#include "lapack_cpp/base/enums.hpp"
#include "lapack_cpp/base/memory_block.hpp"

namespace lapack_cpp {

// Forward declaration of ConstVector
template <typename T, typename idx_t = lapack_idx_t>
class ConstVector;

/**
 * A vector class that can be used to store vectors of arbitrary size.
 *
 * @tparam T this is a template parameter that specifies the type of the
 *           elements of the vector.
 */
template <typename T, typename idx_t = lapack_idx_t>
class Vector {
   public:
    typedef T value_type;
    typedef idx_t index_type;

    // Constructor for a vector of size n
    template <bool aligned>
    Vector(idx_t n,
           const MemoryBlock<T, idx_t, aligned>& m_block,
           idx_t stride = 1)
        : n_(n), stride_(stride), data_(m_block.ptr())
    {
        assert(n >= 0);
        assert(m_block.size() >= n);
    }

    // Constructor for a vector of size n
    Vector(idx_t n, T* data, idx_t stride = 1)
        : n_(n), stride_(stride), data_(data)
    {
        assert(n >= 0);
    }

    // Copy constructor
    // Note: the default copy constructor would probably work, but we want to
    // emphasize that we are doing a shallow copy here.
    Vector(const Vector<T>& v) : n_(v.n_), stride_(v.stride_), data_(v.data_) {}

    // Assignment operator
    void operator=(const Vector<T>& v)
    {
        assert(v.size() == this->size());
        for (int i = 0; i < this->size(); ++i) {
            (*this)[i] = v[i];
        }
    }

    // Make a const promotion of the vector
    inline ConstVector<T, idx_t> as_const() const
    {
        return ConstVector<T, idx_t>(n_, data_, stride_);
    }

    // Returns the size of the vector
    inline idx_t size() const { return n_; }

    inline idx_t stride() const { return stride_; }

    // Returns the i-th element of the vector
    // Note: this function will only be called when the vector is not const
    // therefore, we return the element by reference.
    inline T& operator[](const idx_t i) const
    {
        assert(i < n_);
        return data_[i * stride_];
    }

    /**
     * The subvector starts at index start and ends at index end - 1.
     * The stride is the step size between two elements of the subvector.
     */
    inline Vector subvector(idx_t start, idx_t end, idx_t stride = 1) const
    {
        assert(start >= 0);
        assert(end <= n_);
        assert(start < end);
        assert(stride > 0);
        return Vector((end - start) / stride, &data_[start * stride_],
                      stride * stride_);
    }

    /**
     * Return a pointer to the data of the vector.
     */
    inline T* ptr() const { return data_; }

   private:
    const idx_t n_;
    const idx_t stride_;
    T* data_;
};

/**
 * A vector class that can be used to store vectors of arbitrary size.
 *
 * @tparam T this is a template parameter that specifies the type of the
 *           elements of the vector.
 */
template <typename T, typename idx_t>
class ConstVector {
   public:
    typedef T value_type;
    typedef idx_t index_type;
    // Constructor for a vector of size n
    template <bool aligned>
    ConstVector(idx_t n,
                const MemoryBlock<T, idx_t, aligned>& m_block,
                idx_t stride = 1)
        : n_(n), stride_(stride), data_(m_block.ptr())
    {
        assert(n >= 0);
        assert(m_block.size() >= n);
    }

    // Constructor for a vector of size n
    ConstVector(idx_t n, const T* data, idx_t stride = 1)
        : n_(n), stride_(stride), data_(data)
    {
        assert(n >= 0);
    }

    // Copy constructor
    // Note: the default copy constructor would probably work, but we want to
    // emphasize that we are doing a shallow copy here.
    ConstVector(const ConstVector<T, idx_t>& v)
        : n_(v.n_), stride_(v.stride_), data_(v.data_)
    {}

    // Const promotion constructor
    ConstVector(const Vector<T, idx_t>& v)
        : n_(v.size()), stride_(v.stride()), data_(v.ptr())
    {}

    // Assignment operator
    void operator=(const Vector<T>& v)
    {
        assert(v.size() == this->size());
        for (int i = 0; i < this->size(); ++i) {
            (*this)[i] = v[i];
        }
    }

    // Returns the size of the vector
    inline idx_t size() const { return n_; }

    inline idx_t stride() const { return stride_; }

    // Returns the i-th element of the vector
    // Note: this function will only be called when the vector is not const
    // therefore, we return the element by reference.
    inline T operator[](const idx_t i) const
    {
        assert(i < n_);
        return data_[i * stride_];
    }

    /**
     * The subvector starts at index start and ends at index end - 1.
     * The stride is the step size between two elements of the subvector.
     *
     * Note: this function will also be called if the matrix is const, even
     * though the returned object is not const. Solving this problem would
     * require either an inefficient or a more complex design.
     */
    inline ConstVector subvector(idx_t start, idx_t end, idx_t stride = 1) const
    {
        assert(start >= 0);
        assert(end <= n_);
        assert(start < end);
        assert(stride > 0);
        return ConstVector((end - start) / stride, &data_[start * stride_],
                           stride * stride_);
    }

    /**
     * Return a pointer to the data of the vector.
     */
    inline const T* ptr() const { return data_; }

   private:
    const idx_t n_;
    const idx_t stride_;
    const T* data_;
};

}  // namespace lapack_cpp

#endif