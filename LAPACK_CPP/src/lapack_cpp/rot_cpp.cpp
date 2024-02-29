#ifdef USE_FORTRAN_BLAS
    #include <cassert>
    #include <complex>
    #include <type_traits>

    #include "lapack_c.h"
    #include "lapack_cpp/base.hpp"
    #include "lapack_cpp/blas/rot.hpp"

namespace lapack_cpp {

/**
 * Templated wrapper around c functions, will be instantiated for each type.
 *
 * @tparam T
 * @param A
 * @param B
 * @param C
 */
template <typename T, typename TC, typename TS, typename idx_t>
inline void rot_c_wrapper(const Vector<T, idx_t>& x,
                          const Vector<T, idx_t>& y,
                          TC c,
                          TS s)
{
    assert(x.size() == y.size());

    if constexpr (std::is_same<T, double>::value) {
        lapack_c_drot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
    }
    else if constexpr (std::is_same<T, float>::value) {
        lapack_c_srot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
    }
    else if constexpr (std::is_same<T, std::complex<float>>::value) {
        lapack_c_crot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
    }
    else if constexpr (std::is_same<T, std::complex<double>>::value) {
        lapack_c_zrot(x.size(), x.ptr(), x.stride(), y.ptr(), y.stride(), c, s);
    }
    else {
        assert(false);
    }
}

    // We have a lot of types to instantiate for, so we use a macro to avoid
    // repetition.
    #define INSTANTIATE_ROT(T, TC, idx_t)                                    \
        template <>                                                          \
        void rot(const Vector<T, idx_t>& x, const Vector<T, idx_t>& y, TC c, \
                 T s)                                                        \
        {                                                                    \
            rot_c_wrapper(x, y, c, s);                                       \
        }

INSTANTIATE_ROT(float, float, lapack_idx_t)
INSTANTIATE_ROT(double, double, lapack_idx_t)
INSTANTIATE_ROT(std::complex<float>, float, lapack_idx_t)
INSTANTIATE_ROT(std::complex<double>, double, lapack_idx_t)
INSTANTIATE_ROT(float, float, int)
INSTANTIATE_ROT(double, double, int)
INSTANTIATE_ROT(std::complex<float>, float, int)
INSTANTIATE_ROT(std::complex<double>, double, int)

    #undef INSTANTIATE_ROT

}  // namespace lapack_cpp
#endif  // USE_FORTRAN_BLAS