#ifndef LAPACK_CPP_LARTG_HPP
#define LAPACK_CPP_LARTG_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Generate a plane rotation.
 */
template <typename T>
void lartg(T f, T g, real_t<T>& c, T& s, T& r);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_GEMV_HPP