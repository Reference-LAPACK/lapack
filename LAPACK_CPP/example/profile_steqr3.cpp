
#include <chrono> // for high_resolution_clock
#include <complex>
#include <iostream>

#include "lapack_cpp.hpp"

using namespace lapack_cpp;

template <typename T,
          Layout layout = Layout::ColMajor,
          typename idx_t = lapack_idx_t>
void profile_steqr3(idx_t n, bool use_fortran = false)
{
    MemoryBlock<T, idx_t> Z_(n, n, layout);
    Matrix<T, layout, idx_t> Z(n, n, Z_);

    MemoryBlock<T, idx_t> d_(n);
    Vector<T, idx_t> d(n, d_);

    MemoryBlock<T, idx_t> e_(n - 1);
    Vector<T, idx_t> e(n - 1, e_);

    MemoryBlock<T, idx_t> d_copy_(n);
    Vector<T, idx_t> d_copy(n, d_copy_);

    MemoryBlock<T, idx_t> e_copy_(n - 1);
    Vector<T, idx_t> e_copy(n - 1, e_copy_);

    randomize(d);
    randomize(e);

    d_copy = d;
    e_copy = e;

    CompQ compz = CompQ::Initialize;

    const idx_t n_timings = 100;
    const idx_t n_warmup = 50;
    std::vector<float> timings(n_timings);

    if (use_fortran)
    {
        MemoryBlock<real_t<T>, idx_t> rwork(2 * n - 2);
        for (idx_t i = 0; i < n_timings; ++i)
        {
            d = d_copy;
            e = e_copy;
            auto start = std::chrono::high_resolution_clock::now();
            steqr(compz, d, e, Z, rwork);
            auto end = std::chrono::high_resolution_clock::now();
            timings[i] = std::chrono::duration<float>(end - start).count();
            std::cout << i << " " << timings[i] << '\r' << std::flush;
        }
    }
    else
    {
        MemoryBlock<T, idx_t> work(steqr3_workquery(compz, d, e, Z));
        MemoryBlock<real_t<T>, idx_t> rwork(steqr3_rworkquery(compz, d, e, Z));
        for (idx_t i = 0; i < n_timings; ++i)
        {
            d = d_copy;
            e = e_copy;
            auto start = std::chrono::high_resolution_clock::now();
            steqr3(compz, d, e, Z, work, rwork);
            auto end = std::chrono::high_resolution_clock::now();
            timings[i] = std::chrono::duration<float>(end - start).count();
            std::cout << i << " " << timings[i] << '\r' << std::flush;
        }
    }

    float mean = 0;
    for (idx_t i = n_warmup; i < n_timings; ++i)
        mean += timings[i];
    mean /= n_timings - n_warmup;

    float std_dev = 0;
    for (idx_t i = n_warmup; i < n_timings; ++i)
        std_dev += (timings[i] - mean) * (timings[i] - mean);
    std_dev = std::sqrt(std_dev / (n_timings - n_warmup - 1));

    std::cout << "n = " << n << ", using fortran: " << use_fortran << ", mean time = " << mean
              << " s"
              << ", std dev = " << std_dev / mean * 100 << " %"
              << std::endl;
}

int main()
{
    typedef lapack_idx_t idx_t;
    typedef double T;

    for (int n = 32; n <= 2000; n *= 2)
    {
        profile_steqr3<T, Layout::ColMajor, idx_t>(n, false);
    }
    for (int n = 32; n <= 2000; n *= 2)
    {
        profile_steqr3<T, Layout::ColMajor, idx_t>(n, true);
    }
    return 0;
}