#ifndef LAPACK_CPP_STEQR_HPP
#define LAPACK_CPP_STEQR_HPP

#include "lapack_cpp/base.hpp"

namespace lapack_cpp
{

    /**
     * Computes all eigenvalues and, optionally, eigenvectors of a symmetric tridiagonal matrix using the QR algorithm.
     */
    template <typename T,
              Layout layout,
              typename idx_t,
              bool aligned>
    int steqr(
        CompQ compz,
        Vector<real_t<T>, idx_t> d,
        Vector<real_t<T>, idx_t> e,
        Matrix<T, layout, idx_t> Z,
        MemoryBlock<real_t<T>, idx_t, aligned> &work);

}

#endif // LAPACK_CPP_STEQR_HPP