#ifndef LAPACK_CPP_STEQR3_HPP
#define LAPACK_CPP_STEQR3_HPP

#include "lapack_cpp/base.hpp"
#include "lapack_cpp/lapack/lasr3.hpp"
#include "lapack_cpp/lapack/lartg.hpp"
#include "lapack_cpp/lapack/lasrt.hpp"
#include "lapack_cpp/lapack/laev2.hpp"
#include "lapack_cpp/lapack/lae2.hpp"
#include "lapack_cpp/lapack/lapy2.hpp"


namespace lapack_cpp
{
    template <typename T,
              Layout layout,
              typename idx_t>
    idx_t steqr3_workquery(CompQ compz,
                           Vector<real_t<T>, idx_t> d,
                           Vector<real_t<T>, idx_t> e,
                           Matrix<T, layout, idx_t> Z)
    {
        if( compz == CompQ::No )
            return 0;
        const auto nb = idx_t{32};
        const idx_t n = d.size();

        assert(n == e.size() + 1);
        assert(n == Z.num_columns());
        assert(n == Z.num_rows());

        // Create empty wrappers for C and S to do workspace query
        ConstMatrix<real_t<T>, layout, idx_t> C(n - 1, nb, nullptr);
        ConstMatrix<real_t<T>, layout, idx_t> S(n - 1, nb, nullptr);

        return lasr3_workquery(Side::Right, Direction::Forward, C, S, Z);
    }

    template <typename T,
              Layout layout,
              typename idx_t>
    idx_t steqr3_rworkquery(CompQ compz,
                            Vector<real_t<T>, idx_t> d,
                            Vector<real_t<T>, idx_t> e,
                            Matrix<T, layout, idx_t> Z)
    {
        const auto nb = idx_t{32};
        const idx_t n = d.size();

        assert(n == e.size() + 1);
        assert(n == Z.num_columns());
        assert(n == Z.num_rows());

        if (layout == Layout::ColMajor)
            return 2 * calc_ld<real_t<T>, idx_t>(n - 1) * nb;
        else
            return 2 * (n - 1) * calc_ld<real_t<T>, idx_t>(nb);
    }

    /**
     * STEQR3 computes all eigenvalues and, optionally, eigenvectors of a
     * hermitian tridiagonal matrix using the implicit QL or QR method.
     * The eigenvectors of a full or band hermitian matrix can also be found
     * if SYTRD or SPTRD or SBTRD has been used to reduce this matrix to
     * tridiagonal form.
     *
     * STEQR3 is a variant of STEQR that uses lasr3 to efficiently
     * accumulate the rotations into the eigenvector matrix.
     *
     * @return  0 if success
     *
     * @param[in] want_z bool
     *
     * @param[in,out] d Real vector of length n.
     *      On entry, diagonal elements of the bidiagonal matrix B.
     *      On exit, the singular values of B in decreasing order.
     *
     * @param[in,out] e Real vector of length n-1.
     *      On entry, off-diagonal elements of the bidiagonal matrix B.
     *      On exit, the singular values of B in decreasing order.
     *
     * @param[in,out] Z n-by-m matrix.
     *      On entry, the n-by-n unitary matrix used in the reduction
     *      to tridiagonal form.
     *      On exit, if info = 0, and want_z=true then Z contains the
     *      orthonormal eigenvectors of the original Hermitian matrix.
     * 
     * @param[out] work Workspace array whose size is returned by
     *    steqr3_workquery.
     * 
     * @param[out] rwork Real workspace array whose size is returned by
     *     steqr3_rworkquery.
     *
     * @ingroup computational
     */
    template <typename T,
              Layout layout,
              typename idx_t,
              bool aligned>
    int steqr3(
        CompQ compz,
        Vector<real_t<T>, idx_t> d,
        Vector<real_t<T>, idx_t> e,
        Matrix<T, layout, idx_t> Z,
        MemoryBlock<T, idx_t, aligned> &work,
        MemoryBlock<real_t<T>, idx_t, aligned> &rwork)
    {
        using T_real = real_t<T>;

        // constants
        const T_real two(2);
        const T_real one(1);
        const T_real zero(0);
        const idx_t n = d.size();

        assert(n == e.size() + 1);
        assert(n == Z.num_columns());
        assert(n == Z.num_rows());
        assert( rwork.size() >= steqr3_rworkquery(want_z, d, e, Z) );
        assert( work.size() >= steqr3_workquery(want_z, d, e, Z) );

        // Amount of rotation sequences to generate before applying it.
        const auto nb = idx_t{32};

        // Quick return if possible
        if (n == 0)
            return 0;
        if (n == 1)
        {
            if (compz == CompQ::Initialize)
                Z(0, 0) = one;
            return 0;
        }

        bool want_z = compz != CompQ::No;

        if( compz == CompQ::Initialize )
        {
            for (idx_t j = 0; j < n; ++j)
            {
                for (idx_t i = 0; i < n; ++i)
                {
                    Z(i, j) = zero;
                }
                Z(j, j) = one;
            }
        }

        // Determine the unit roundoff and over/underflow thresholds.
        const T_real eps = ulp<T_real>();
        const T_real eps2 = eps * eps;
        const T_real safmin = safe_min<T_real>();

        // Compute the eigenvalues and eigenvectors of the tridiagonal
        // matrix.
        const idx_t itmax = 30 * n;

        // istart and istop determine the active block
        idx_t istart = 0;
        idx_t istop = n;

        // Keep track of previous istart and istop to know when to change direction
        idx_t istart_old = -1;
        idx_t istop_old = -1;

        // Keep track of the delayed rotation sequences
        idx_t i_block = 0;

        // C and S are used to store the rotation sequences
        Matrix<T_real, layout, idx_t> C(n - 1, nb, rwork);
        auto rwork2 = rwork.remainder(n - 1, nb);
        Matrix<T_real, layout, idx_t> S(n - 1, nb, rwork2);

        // Initialize C and S to identity rotations
        for (idx_t j = 0; j < nb; ++j)
        {
            for (idx_t i = 0; i < n - 1; ++i)
            {
                C(i, j) = one;
                S(i, j) = zero;
            }
        }

        // Direction to chase the bulge
        // This variable is reevaluated for every new subblock
        Direction direction = d[0] > d[n - 1] ? Direction::Forward : Direction::Backward;

        // Main loop
        for (idx_t iter = 0; iter < itmax; iter++)
        {
            // If we have saved up enough rotations, apply them
            if (want_z and (i_block >= nb or iter == itmax or istop <= 1))
            {
                idx_t i_block2 = std::min<idx_t>(i_block + 1, nb);

                // Find smallest index where rotation is not identity
                idx_t i_start_block = n - 1;
                for (idx_t i = 0; i < i_block2; ++i)
                {
                    for (idx_t j = 0; j < i_start_block; ++j)
                        if (C(j, i) != one or S(j, i) != zero)
                        {
                            i_start_block = j;
                            break;
                        }
                }

                // Find largest index where rotation is not identity
                idx_t i_stop_block = 0;
                for (idx_t i = 0; i < i_block2; ++i)
                {
                    for (idx_t j2 = n - 1; j2 > i_stop_block; --j2)
                    {
                        idx_t j = j2 - 1;
                        if (C(j, i) != one or S(j, i) != zero)
                        {
                            i_stop_block = j;
                            break;
                        }
                    }
                }

                if (i_start_block < i_stop_block + 1)
                {
                    auto C2 = C.submatrix(i_start_block, i_stop_block + 1, 0, i_block2);
                    auto S2 = S.submatrix(i_start_block, i_stop_block + 1, 0, i_block2);
                    auto Z2 = Z.submatrix(0, n, i_start_block, i_stop_block + 2);
                    lasr3(Side::Right, direction, C2.as_const(), S2.as_const(), Z2, work);
                }
                // Reset block
                i_block = 0;

                // Initialize C and S to identity rotations
                for (idx_t j = 0; j < nb; ++j)
                {
                    for (idx_t i = 0; i < n - 1; ++i)
                    {
                        C(i, j) = one;
                        S(i, j) = zero;
                    }
                }
            }
            if (iter == itmax)
            {
                // The QR algorithm failed to converge, return with error.
                return istop;
            }

            if (istop <= 1)
            {
                // All eigenvalues have been found, exit and return 0.
                break;
            }

            // Find active block
            for (idx_t i = istop - 1; i > istart; --i)
            {
                if ((e[i - 1] * e[i - 1]) <=
                    (eps2 * abs(d[i - 1])) * abs(d[i]) + safmin)
                {
                    e[i - 1] = zero;
                    istart = i;
                    break;
                }
            }

            // An eigenvalue has split off, reduce istop and start the loop again
            if (istart == istop - 1)
            {
                istop = istop - 1;
                istart = 0;
                continue;
            }

            // A 2x2 block has split off, handle separately
            if (istart + 1 == istop - 1)
            {
                T_real s1, s2;
                if (want_z)
                {
                    T_real cs, sn;
                    laev2(d[istart], e[istart], d[istart + 1], s1, s2, cs, sn);

                    // Store rotation
                    C(istart, i_block) = cs;
                    S(istart, i_block) = sn;
                    // Normally, we would increment i_block here, but we do not
                    // because the block has been deflated and should not interfere
                    // with the next sequence.
                }
                else
                {
                    lae2(d[istart], e[istart], d[istart + 1], s1, s2);
                }
                d[istart] = s1;
                d[istart + 1] = s2;
                e[istart] = zero;

                istop = istop - 2;
                istart = 0;
                continue;
            }

            // Choose between QL and QR iteration
            if (istart >= istop_old or istop <= istart_old)
            {
                Direction direction_new;
                // We prefer to keep the current direction, because switching
                // direction forces us to apply the rotations in the block, which
                // may be inefficient if the block is small.
                if (direction == Direction::Forward)
                {
                    direction_new = 10. * abs(d[istart]) > abs(d[istop - 1]) ? Direction::Forward : Direction::Backward;
                }
                else
                {
                    direction_new =
                        abs(d[istart]) > 10. * abs(d[istop - 1]) ? Direction::Forward : Direction::Backward;
                }

                if (direction_new != direction)
                {
                    // We don't want different directions in our saved rotation matrices
                    // So we apply all the current rotations and reset the block
                    idx_t i_block2 = std::min<idx_t>(i_block + 1, nb);

                    // Find smallest index where rotation is not identity
                    idx_t i_start_block = n - 1;
                    for (idx_t i = 0; i < i_block2; ++i)
                    {
                        for (idx_t j = 0; j < i_start_block; ++j)
                            if (C(j, i) != one or S(j, i) != zero)
                            {
                                i_start_block = j;
                                break;
                            }
                    }

                    // Find largest index where rotation is not identity
                    idx_t i_stop_block = 0;
                    for (idx_t i = 0; i < i_block2; ++i)
                    {
                        for (idx_t j2 = n - 1; j2 > i_stop_block; --j2)
                        {
                            idx_t j = j2 - 1;
                            if (C(j, i) != one or S(j, i) != zero)
                            {
                                i_stop_block = j;
                                break;
                            }
                        }
                    }

                    if (i_start_block < i_stop_block + 1)
                    {
                        auto C2 = C.submatrix(i_start_block, i_stop_block + 1, 0, i_block2);
                        auto S2 = S.submatrix(i_start_block, i_stop_block + 1, 0, i_block2);
                        auto Z2 = Z.submatrix(0, n, i_start_block, i_stop_block + 2);

                        lasr3(Side::Right, direction, C2.as_const(), S2.as_const(), Z2, work);
                    }
                    // Reset block
                    i_block = 0;

                    // Initialize C and S to identity rotations
                    for (idx_t j = 0; j < nb; ++j)
                    {
                        for (idx_t i = 0; i < n - 1; ++i)
                        {
                            C(i, j) = one;
                            S(i, j) = zero;
                        }
                    }
                }
                direction = direction_new;
            }
            istart_old = istart;
            istop_old = istop;

            if (direction == Direction::Forward)
            {
                // QR iteration

                // Form shift using last 2x2 block of the active matrix
                T_real p = d[istop - 1];
                T_real g = (d[istop - 2] - p) / (two * e[istop - 2]);
                T_real r = lapy2(g, one);
                g = d[istart] - p + e[istop - 2] / (T_real)(g + (sgn(g) * r));

                T_real s = one;
                T_real c = one;
                p = zero;

                // Chase bulge from top to bottom
                for (idx_t i = istart; i < istop - 1; ++i)
                {
                    T_real f = s * e[i];
                    T_real b = c * e[i];
                    lartg(g, f, c, s, r);
                    if (i != istart)
                        e[i - 1] = r;
                    g = d[i] - p;
                    r = (d[i + 1] - g) * s + two * c * b;
                    p = s * r;
                    d[i] = g + p;
                    g = c * r - b;
                    // If eigenvalues are desired, then apply rotations
                    if (want_z)
                    {
                        // Store rotation for later
                        C(i, i_block) = c;
                        S(i, i_block) = s;
                    }
                }
                d[istop - 1] = d[istop - 1] - p;
                e[istop - 2] = g;
            }
            else
            {
                // QL iteration

                // Form shift using last 2x2 block of the active matrix
                T_real p = d[istart];
                T_real g = (d[istart + 1] - p) / (two * e[istart]);
                T_real r = lapy2(g, one);
                g = d[istop - 1] - p + e[istart] / (T_real)(g + (sgn(g) * r));

                T_real s = one;
                T_real c = one;
                p = zero;

                // Chase bulge from bottom to top
                for (idx_t i = istop - 1; i > istart; --i)
                {
                    T_real f = s * e[i - 1];
                    T_real b = c * e[i - 1];
                    lartg(g, f, c, s, r);
                    if (i != istop - 1)
                        e[i] = r;
                    g = d[i] - p;
                    r = (d[i - 1] - g) * s + two * c * b;
                    p = s * r;
                    d[i] = g + p;
                    g = c * r - b;
                    // If eigenvalues are desired, then apply rotations
                    if (want_z)
                    {
                        // Store rotation for later
                        C(i - 1, i_block) = c;
                        S(i - 1, i_block) = -s;
                    }
                }
                d[istart] = d[istart] - p;
                e[istart] = g;
            }
            i_block++;
        }

        // Order eigenvalues and eigenvectors
        if (!want_z)
        {
            // Order eigenvalues
            lasrt(IncDec::Increasing, d);
        }
        else
        {
            // Use selection sort to minize swaps of eigenvectors
            for (idx_t i = 0; i < n - 1; ++i)
            {
                idx_t k = i;
                T_real p = d[i];
                for (idx_t j = i + 1; j < n; ++j)
                {
                    if (d[j] < p)
                    {
                        k = j;
                        p = d[j];
                    }
                }
                if (k != i)
                {
                    d[k] = d[i];
                    d[i] = p;
                    auto z1 = Z.column(i);
                    auto z2 = Z.column(k);
                    for (idx_t j = 0; j < n; ++j)
                    {
                        auto temp = z1[j];
                        z1[j] = z2[j];
                        z2[j] = temp;
                    }
                }
            }
        }

        return 0;
    }

}

#endif // LAPACK_CPP_STEQR3_HPP