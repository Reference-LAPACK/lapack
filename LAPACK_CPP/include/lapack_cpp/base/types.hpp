#ifndef LAPACK_CPP_BASE_TYPES_HPP
#define LAPACK_CPP_BASE_TYPES_HPP

#include <complex>
#include <cstdint>
#include <type_traits>

// While the library is templated and can handle different integers,
// this is the default integer type used by the library.
typedef int64_t lapack_idx_t;

// declare conj for real types
template <typename T,
          std::enable_if_t<std::is_floating_point_v<T>, bool> = true>
inline T conj(const T& x)
{
    return x;
}

template <typename T>
struct real_type_struct {
    typedef T type;
};

template <typename T>
struct real_type_struct<std::complex<T>> {
    typedef T type;
};

template <typename T>
using real_t = typename real_type_struct<T>::type;


template <typename T>
constexpr real_t<T> ulp() noexcept
{
    return std::numeric_limits<real_t<T>>::epsilon();
}

template <typename T>
constexpr real_t<T> safe_min() noexcept
{
    const int fradix = std::numeric_limits<real_t<T>>::radix;
    const int expm = std::numeric_limits<real_t<T>>::min_exponent;
    const int expM = std::numeric_limits<real_t<T>>::max_exponent;

    return pow(fradix, real_t<T>(std::max(expm - 1, 1 - expM)));
}

/// Sign function for real numbers
/// Note, this has the following behavior:
/// sgn(x) = 1 if x > 0
/// sgn(x) = -1 if x < 0
/// sgn(0) = 1
/// sgn(-0) = 1
/// sgn(+Inf) = 1
/// sgn(-Inf) = 1
/// sgn(NaN) = -1
template <typename T>
constexpr T sgn(const T& val)
{
    return (val >= T(0)) ? T(1) : T(-1);
}

#endif  // LAPACK_CPP_BASE_TYPES_HPP