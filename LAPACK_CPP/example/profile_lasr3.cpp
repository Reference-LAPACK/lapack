
#include <chrono>  // for high_resolution_clock
#include <complex>
#include <iostream>

#include "lapack_cpp.hpp"

using namespace lapack_cpp;

template <typename T,
          Layout layout = Layout::ColMajor,
          typename idx_t = lapack_idx_t>
void profile_lasr3(idx_t m, idx_t n, idx_t k, Side side, Direction direction)
{
    idx_t n_rot = side == Side::Left ? m - 1 : n - 1;

    MemoryBlock<T, idx_t> A_(m, n, layout);
    Matrix<T, layout, idx_t> A(m, n, A_);

    MemoryBlock<real_t<T>, idx_t> C_(n_rot, k, layout);
    Matrix<real_t<T>, layout, idx_t> C(n_rot, k, C_);

    MemoryBlock<T, idx_t> S_(n_rot, k, layout);
    Matrix<T, layout, idx_t> S(n_rot, k, S_);

    randomize(A);

    // Generate random, but valid rotations
    randomize(C);
    randomize(S);
    for (idx_t i = 0; i < n_rot; ++i) {
        for (idx_t j = 0; j < k; ++j) {
            T f = C(i, j);
            T g = S(i, j);
            T r;
            lartg(f, g, C(i, j), S(i, j), r);
        }
    }

    MemoryBlock<T, idx_t> A_copy_(m, n, layout);
    Matrix<T, layout, idx_t> A_copy(m, n, A_copy_);

    A_copy = A;

    MemoryBlock<T, idx_t> work(
        lasr3_workquery(side, direction, C.as_const(), S.as_const(), A));

    const idx_t n_timings = 100;
    const idx_t n_warmup = 20;
    std::vector<float> timings(n_timings);

    for (idx_t i = 0; i < n_timings; ++i) {
        A = A_copy;
        auto start = std::chrono::high_resolution_clock::now();
        lasr3(side, direction, C.as_const(), S.as_const(), A, work);
        auto end = std::chrono::high_resolution_clock::now();
        timings[i] = std::chrono::duration<float>(end - start).count();
    }

    float mean = 0;
    for (idx_t i = n_warmup; i < n_timings; ++i)
        mean += timings[i];
    mean /= n_timings - n_warmup;

    long nflops =
        side == Side::Left ? 6. * (m - 1) * n * k : 6. * (n - 1) * m * k;

    std::cout << "m = " << m << ", n = " << n << ", k = " << k
              << ", side = " << (char)side
              << ", direction = " << (char)direction << ", mean time = " << mean
              << " s"
              << ", flop rate = " << nflops / mean * 1.0e-9 << " GFlops"
              << std::endl;
}

int main()
{
    typedef lapack_idx_t idx_t;
    typedef double T;

    const idx_t nb = 1000;
    const idx_t k = 64;

    for (int s = 0; s < 2; ++s) {
        Side side = s == 0 ? Side::Left : Side::Right;
        for (int d = 0; d < 2; ++d) {
            Direction direction =
                d == 0 ? Direction::Forward : Direction::Backward;
            for (idx_t n_rot = 99; n_rot <= 1600; n_rot += 100) {
                idx_t m = side == Side::Left ? n_rot + 1 : nb;
                idx_t n = side == Side::Left ? nb : n_rot + 1;

                profile_lasr3<T, Layout::ColMajor, lapack_idx_t>(m, n, k, side,
                                                                 direction);
            }
        }
    }
    return 0;
}