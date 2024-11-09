#ifndef LAPACK_CPP_LAEV2_HPP
#define LAPACK_CPP_LAEV2_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"

namespace lapack_cpp {

/** Computes the eigenvalues and eigenvector of a real symmetric 2x2 matrix A
 *  [ a b ]
 *  [ b c ]
 *  On exit, the decomposition satisfies:
 *  [ cs  sn ] [ a b ] [ cs -sn ] = [ s1  0  ]
 *  [ -sn cs ] [ b c ] [ sn  cs ]   [ 0   s2 ]
 *  where cs*cs + sn*sn = 1.
 *
 * @param[in] a
 *      Element (0,0) of A.
 * @param[in] b
 *      Element (0,1) and (1,0) of A.
 * @param[in] c
 *      Element (1,1) of A.
 * @param[out] s1
 *      The eigenvalue of A with the largest absolute value.
 * @param[out] s2
 *      The eigenvalue of A with the smallest absolute value.
 * @param[out] cs
 *      The cosine of the rotation matrix.
 * @param[out] sn
 *      The sine of the rotation matrix.
 *
 * \verbatim
 *  s1 is accurate to a few ulps barring over/underflow.
 *
 *  s2 may be inaccurate if there is massive cancellation in the
 *  determinant a*c-b*b; higher precision or correctly rounded or
 *  correctly truncated arithmetic would be needed to compute s2
 *  accurately in all cases.
 *
 *  cs and sn are accurate to a few ulps barring over/underflow.
 *
 *  Overflow is possible only if s1 is within a factor of 5 of overflow.
 *  Underflow is harmless if the input data is 0 or exceeds
 *     underflow_threshold / macheps.
 * \endverbatim
 *
 *
 * @ingroup auxiliary
 */
template <typename T>
void laev2(T a, T b, T c, T& s1, T& s2, T& cs, T& sn);

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_LAEV2_HPP