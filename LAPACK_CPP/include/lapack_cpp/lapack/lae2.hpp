#ifndef LAPACK_CPP_LAE2_HPP
#define LAPACK_CPP_LAE2_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"

namespace lapack_cpp {

/** Computes the eigenvalues of a real symmetric 2x2 matrix A
 *  [ a b ]
 *  [ b c ]
 */
template <typename T>
void lae2(T a, T b, T c, T& s1, T& s2);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_LAE2_HPP