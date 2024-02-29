#ifndef LAPACK_CPP_ROT_HPP
#define LAPACK_CPP_ROT_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Apply a plane rotation to a pair of vectors.
 */
template <typename T, typename TS, typename idx_t>
void rot(const Vector<T, idx_t>& x, const Vector<T, idx_t>& y, real_t<T> c, TS s){
    assert(x.size() == y.size());
    const idx_t n = x.size();
    for (idx_t i = 0; i < n; ++i) {
        T temp = -conj(s) * x[i] + c * y[i];
        x[i] = c * x[i] + s * y[i];
        y[i] = temp;
    }
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_GEMV_HPP