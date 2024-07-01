#include <cassert>
#include <complex>
#include <type_traits>

#include "lapack_c.h"
#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/steqr.hpp"

namespace lapack_cpp
{

    template <>
    int steqr(CompQ compz, Vector<float, int> d, Vector<float, int> e, Matrix<float, Layout::ColMajor, int> Z, MemoryBlock<float, int, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_ssteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }

    template <>
    int steqr(CompQ compz, Vector<double, int> d, Vector<double, int> e, Matrix<double, Layout::ColMajor, int> Z, MemoryBlock<double, int, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_dsteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }


    template <>
    int steqr(CompQ compz, Vector<float, int> d, Vector<float, int> e, Matrix<std::complex<float>, Layout::ColMajor, int> Z, MemoryBlock<float, int, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_csteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }

    template <>
    int steqr(CompQ compz, Vector<double, int> d, Vector<double, int> e, Matrix<std::complex<double>, Layout::ColMajor, int> Z, MemoryBlock<double, int, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_zsteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }


    template <>
    int steqr(CompQ compz, Vector<float, int64_t> d, Vector<float, int64_t> e, Matrix<float, Layout::ColMajor, int64_t> Z, MemoryBlock<float, int64_t, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_ssteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }

    template <>
    int steqr(CompQ compz, Vector<double, int64_t> d, Vector<double, int64_t> e, Matrix<double, Layout::ColMajor, int64_t> Z, MemoryBlock<double, int64_t, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_dsteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }


    template <>
    int steqr(CompQ compz, Vector<float, int64_t> d, Vector<float, int64_t> e, Matrix<std::complex<float>, Layout::ColMajor, int64_t> Z, MemoryBlock<float, int64_t, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_csteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }

    template <>
    int steqr(CompQ compz, Vector<double, int64_t> d, Vector<double, int64_t> e, Matrix<std::complex<double>, Layout::ColMajor, int64_t> Z, MemoryBlock<double, int64_t, true> &work)
    {
        assert(d.size() >= 0);
        assert(d.size() == e.size() + 1);
        assert(d.size() == Z.num_rows());
        assert(d.size() == Z.num_columns());
        assert(2*d.size()-2 == work.size());
        // fortran code only supports stride 1
        assert(d.stride() == 1);
        assert(e.stride() == 1);
        lapack_idx info;
        lapack_c_zsteqr((char)compz, d.size(), d.ptr(), e.ptr(), Z.ptr(), Z.ldim(), work.ptr(), &info);

        return info;
    }

} // namespace lapack_cpp