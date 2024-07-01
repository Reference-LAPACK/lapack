#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/lasrt.hpp"

namespace lapack_cpp {

template <>
void lasrt<float, int>(IncDec id, Vector<float, int> d)
{
    assert(d.size() >= 0);
    // fortran code only supports stride 1
    assert(d.stride() == 1);
    lapack_idx info;
    lapack_c_slasrt((char) id, d.size(), d.ptr(), &info);
}

template <>
void lasrt<double, int>(IncDec id, Vector<double, int> d)
{
    assert(d.size() >= 0);
    // fortran code only supports stride 1
    assert(d.stride() == 1);
    lapack_idx info;
    lapack_c_dlasrt((char) id, d.size(), d.ptr(), &info);
}

template <>
void lasrt<float, int64_t>(IncDec id, Vector<float, int64_t> d)
{
    assert(d.size() >= 0);
    // fortran code only supports stride 1
    assert(d.stride() == 1);
    lapack_idx info;
    lapack_c_slasrt((char) id, d.size(), d.ptr(), &info);
}

template <>
void lasrt<double, int64_t>(IncDec id, Vector<double, int64_t> d)
{
    assert(d.size() >= 0);
    // fortran code only supports stride 1
    assert(d.stride() == 1);
    lapack_idx info;
    lapack_c_dlasrt((char) id, d.size(), d.ptr(), &info);
}

}  // namespace lapack_cpp