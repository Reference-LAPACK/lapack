#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/lartg.hpp"

namespace lapack_cpp {

template <>
void lartg<float>(float f, float g, float& c, float& s, float& r)
{
    lapack_c_slartg(f, g, &c, &s, &r);
}

template <>
void lartg<double>(double f, double g, double& c, double& s, double& r)
{
    lapack_c_dlartg(f, g, &c, &s, &r);
}

template <>
void lartg<std::complex<float>>(std::complex<float> f,
                                std::complex<float> g,
                                float& c,
                                std::complex<float>& s,
                                std::complex<float>& r)
{
    lapack_c_clartg(f, g, &c, &s, &r);
}

template <>
void lartg<std::complex<double>>(std::complex<double> f,
                                 std::complex<double> g,
                                 double& c,
                                 std::complex<double>& s,
                                 std::complex<double>& r)
{
    lapack_c_zlartg(f, g, &c, &s, &r);
}

}  // namespace lapack_cpp