#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/lae2.hpp"

namespace lapack_cpp {

template <>
void lae2<float>(float a, float b, float c, float& s1, float& s2)
{
    lapack_c_slae2(a, b, c, &s1, &s2);
}

template <>
void lae2<double>(double a, double b, double c, double& s1, double& s2)
{
    lapack_c_dlae2(a, b, c, &s1, &s2);
}

}  // namespace lapack_cpp