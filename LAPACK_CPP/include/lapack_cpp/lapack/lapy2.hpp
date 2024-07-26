#ifndef LAPACK_CPP_LAPY2_HPP
#define LAPACK_CPP_LAPY2_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/utils.hpp"
namespace lapack_cpp {

/**
 * Return \f$ \sqrt{f^2 + g^2} \f$ taking care to avoid overflow and underflow.
 */
template <typename T>
T lapy2(T x, T y){
    // constants
    const T one(1);
    const T zero(0);
    const T xabs = abs(x);
    const T yabs = abs(y);

    T w, z;
    if (xabs > yabs) {
        w = xabs;
        z = yabs;
    }
    else {
        w = yabs;
        z = xabs;
    }

    return (z == zero) ? w : w * sqrt(one + (z / w) * (z / w));
}

}  // namespace lapack_cpp

#endif  // LAPACK_CPP_LAPY2_HPP