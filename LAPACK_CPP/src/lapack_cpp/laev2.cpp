#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/laev2.hpp"

namespace lapack_cpp {

template <>
void laev2<float>(float a, float b, float c, float& s1, float& s2, float& cs, float& sn)
{
    lapack_c_slaev2(a, b, c, &s1, &s2, &cs, &sn);
}

template <>
void laev2<double>(double a, double b, double c, double& s1, double& s2, double& cs, double& sn)
{
    lapack_c_dlaev2(a, b, c, &s1, &s2, &cs, &sn);
}

}  // namespace lapack_cpp