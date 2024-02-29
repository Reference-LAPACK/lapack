
#ifndef LAPACK_CPP_MEMORY_BLOCK_HPP
#define LAPACK_CPP_MEMORY_BLOCK_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <type_traits>

#include "lapack_cpp/base/types.hpp"

namespace lapack_cpp {

/**
 * Calculate the leading dimension of a matrix, given the number of rows.
 * Normally, the leading dimension can just be equal to to the number of rows,
 * but for performance reasons, we want to make sure that the leading dimension
 * is a multiple of 64 bytes. That way, if the first column is aligned, they all
 * are.
 *
 * @tparam T
 * @param m
 * @return int
 */
template <typename T, typename idx_t, bool aligned = true>
inline constexpr idx_t calc_ld(const idx_t m)
{
    if (aligned)
        return ((m - 1 + (64 / sizeof(T))) / (64 / sizeof(T))) *
               (64 / sizeof(T));
    else
        return m;
}

/**
 * A class representing a contiguous block of memory of fixed size.
 *
 * This class owns the data. It is responsible for allocating and deallocating
 * the memory.
 *
 * A special note about const. Declaring a Memory block as const does not make
 * the data const. For example:
 * const MemoryBlock<double> A(10);
 * double* data = A.ptr();
 * data[0] = 1.0; // This is allowed
 */
template <typename T, typename idx_t = lapack_idx_t, bool aligned = true>
class MemoryBlock {
   public:
    // Constructor for a block of size n
    MemoryBlock(idx_t n)
        : n_(n),
          data_(aligned ? (T*)aligned_alloc(64, n_ * sizeof(T))
                        : (T*)malloc(n_ * sizeof(T)))
    {
        assert(n >= 0);
        #ifndef NDEBUG
            // When debugging, fill the memory with NaNs
            // This makes it easier to spot uninitialized memory
            std::fill(data_, data_ + n_, NAN);
        #endif
    }

    // Constructor for a block of sufficient size to store a matrix of size m x
    // n
    MemoryBlock(idx_t m, idx_t n, Layout layout = Layout::ColMajor)
        : n_(layout == Layout::ColMajor ? (calc_ld<T, idx_t, aligned>(m) * n) : (calc_ld<T, idx_t, aligned>(n) * m)),
          data_(aligned ? (T*)aligned_alloc(64, n_ * sizeof(T))
                        : (T*)malloc(n_ * sizeof(T)))
    {
        assert(m >= 0);
        assert(n >= 0);
        #ifndef NDEBUG
            // When debugging, fill the memory with NaNs
            // This makes it easier to spot uninitialized memory
            std::fill(data_, data_ + n_, NAN);
        #endif
    }

    // Non-owning constructor
    MemoryBlock(idx_t n, T* data) : n_(n), data_(data), owns_data_(false)
    {
        assert(n >= 0);
    }

    // Destructor
    ~MemoryBlock()
    {
        if (owns_data_) free((std::remove_const_t<T>*)data_);
    }

    // Copy constructor
    MemoryBlock(const MemoryBlock& other)
        : n_(other.n_), data_(other.data_), owns_data_(false)
    {}

    // Assignment operator
    void operator=(const MemoryBlock& other)
    {
        assert(other.size() == this->size());
        std::copy(other.data_, other.data_ + n_, data_);
    }

    // Memory block that remains after creating a vector of size n
    MemoryBlock remainder(idx_t n) const
    {
        idx_t size = calc_ld<T, idx_t, aligned>(n);
        assert(size <= n_);
        return MemoryBlock(n_ - n, data_ + n);
    }

    // Memory block that remains after creating a vector of size n
    MemoryBlock remainder(idx_t m, idx_t n) const
    {
        idx_t size = calc_ld<T, idx_t, aligned>(calc_ld<T, idx_t, aligned>(m) * n);
        assert(size <= n_);
        return MemoryBlock(n_ - size, data_ + size);
    }

    // Returns the size of the memory block
    inline idx_t size() const { return n_; }

    /**
     * Return a pointer to the data of the vector.
     */
    inline T* ptr() const { return data_; }

   private:
    const idx_t n_;
    T* data_;
    const bool owns_data_ = true;
};

}  // namespace lapack_cpp

#endif