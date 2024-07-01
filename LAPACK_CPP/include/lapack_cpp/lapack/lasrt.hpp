#ifndef LAPACK_CPP_LASRT_HPP
#define LAPACK_CPP_LASRT_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"

namespace lapack_cpp {

/**
 * Sort the numbers in an array.
 */
template <typename T, typename idx_t>
void lasrt(IncDec id, Vector<T, idx_t> d);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_LASRT_HPP