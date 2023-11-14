*> \brief \b DLAQP2RK computes truncated QR factorization with column pivoting of the matrix block using Level 2 BLAS and overwrites m-by-nrhs matrix B with Q**T * B.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAQP2RK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqp2rk.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqp2rk.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqp2rk.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*      SUBROUTINE DLAQP2RK( M, N, NRHS, IOFFSET, KMAX, ABSTOL, RELTOL,
*     $                     KP1, MAXC2NRM, A, LDA, KF, MAXC2NRMK,
*     $                     RELMAXC2NRMK, JPIV, TAU, VN1, VN2, WORK )
*      IMPLICIT NONE
*
*     .. Scalar Arguments ..
*      INTEGER            IOFFSET, KP1, KF, KMAX, LDA, M, N, NRHS
*      DOUBLE PRECISION   ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK,
*     $                   RELTOL
*     ..
*     .. Array Arguments ..
*      INTEGER            JPIV( * )
*      DOUBLE PRECISION   A( LDA, * ), TAU( * ), VN1( * ), VN2( * ),
*     $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAQP2RK computes a truncated (rank K) or full rank Householder QR
*> factorization with column pivoting of the block A(IOFFSET+1:M,1:N).
*> The routine is calling Level 2 BLAS.  The block A(1:IOFFSET,1:N)
*> is accordingly pivoted, but not factorized.  The routine also
*> overwrites the matrix B block stored in A(IOFFSET+1:M,N+1:N+NRHS)
*> with Q(K)**T * B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of
*>          columns of the matrix B. NRHS >= 0.
*> \endverbatim
*>
*> \param[in] IOFFSET
*> \verbatim
*>          IOFFSET is INTEGER
*>          The number of rows of the matrix A that must be pivoted
*>          but not factorized. IOFFSET also represents the number of
*>          columns of the original matrix that have been factorized
*>          in the previous steps.
*>          IOFFSET >= 0.
*> \endverbatim
*>
*> \param[in] MAXK
*> \verbatim
*>          MAXK is INTEGER
*>
*>          The first factorization stopping criterion.
*>
*>          The maximum number of columns of the matrix A to factorize,
*>          i.e. the maximum factorization rank. MAXK >= 0.
*>
*>          a) If MAXK >= min(M-IOFFSET,N), then this stopping
*>                criterion is not used, factorize columns
*>                depending on ABSTOL and RELTOL.
*>
*>          b) If MAXK = 0, then this stopping criterion is
*>                satisfied on input and the routine exits immediately.
*>                This means that the factorization is not performed,
*>                the matrices A and B are not modified, and
*>                the matrix A is itself the residual.
*> \endverbatim
*>
*> \param[in] ABSTOL
*> \verbatim
*>          ABSTOL is DOUBLE PRECISION, cannot be NaN.
*>
*>          The second factorization stopping criterion.
*>
*>          The absolute tolerance (stopping threshold) for
*>          maximum column 2-norm of the residual matrix R22(K).
*>          The algorithm converges (stops the factorization) when
*>          the maximum column 2-norm of the residual matrix R22(K)
*>          is less than or equal to ABSTOL.
*> \endverbatim
*>
*> \param[in] RELTOL
*> \verbatim
*>          RELTOL is DOUBLE PRECISION, cannot be NaN.
*>
*>          The third factorization stopping criterion.
*>
*>          The tolerance (stopping threshold) for the ratio
*>          abs(R(K+1,K+1))/abs(R(1,1)) of the maximum column 2-norm of
*>          the residual matrix R22(K) and the maximum column 2-norm of
*>          the original matrix A. The algorithm converges (stops the
*>          factorization), when abs(R(K+1,K+1))/abs(R(1,1)) A is less
*>          than or equal to RELTOL.
*>
*>          Here, abs(R(1,1)) is the maximum column 2-norm of the
*>          original matrix A; EPS = DLAMCH('E').
*> \endverbatim
*>
*> \param[in] KP1
*> \verbatim
*>          KP1 is INTEGER
*>          The index of the column with the maximum column 2-norm in
*>          the whole original matrix A. KP1 > 0.
*> \endverbatim
*>
*> \param[in] MAXC2NRM
*> \verbatim
*>          MAXC2NRM is DOUBLE PRECISION
*>          The maximum column 2-norm of the whole original matrix.
*>          MAXC2NRMK >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N+NRHS)
*>          On entry:
*>              the M-by-N matrix A and M-by-NRHS matrix B, as in
*>
*>                                  N     NRHS
*>              array_A   =   M  [ mat_A, mat_B ]
*>
*>          On exit:
*>          1. The elements in block A(IOFFSET+1:M,1:KF) below
*>             the diagonal,together with the array TAU, represent
*>             the orthogonal matrix Q(K) as a product of elementary
*>             reflectors.
*>          2. The block of the matrix A stored in A(IOFFSET+1:M,1:KF)
*>             is the triangular factor obtained.
*>          3. The block of the the matrix A stored in A(1:IOFFSET,1:N)
*>             has been accordingly pivoted, but no factorized.
*>          4. The rest of the array A, block A(IOFFSET+1:M,KF+1:N+NRHS).
*>             The left part A(IOFFSET+1:M,KF+1:N) of
*>             this block contains the residual of the matrix A, and
*>             the right part of the block A(IOFFSET+1:M,N+1:N+NRHS)
*>             contains the block of the right-hand-side matrix B. Both
*>             these blocks have been updated by multiplication from
*>             the left by Q**T.
*> \endverbatim
*>
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] KF
*> \verbatim
*>          KF is INTEGER
*>          The number of columns actually factorized.
*> \endverbatim
*>
*> \param[out] MAXC2NRMK
*> \verbatim
*>          MAXC2NRMK is DOUBLE PRECISION
*>          The maximum column 2-norm of the residual matrix A22,
*>          when factorization stopped. MAXC2NRMK >= 0.
*> \endverbatim
*>
*> \param[out] RELMAXC2NRMK
*> \verbatim
*>          RELMAXC2NRMK is DOUBLE PRECISION
*>          The ratio MAXC2NRMK / MAXC2NRM of the maximum column
*>          2-norm of the residual matrix A22 ( when factorization
*>          stopped) and the maximum column 2-norm of the
*>          original matrix A. RELMAXC2NRMK >= 0.
*> \endverbatim
*>
*> \param[out] JPIV
*> \verbatim
*>          JPIV is INTEGER array, dimension (N)
*>          Column pivot indices, for 1 <= j <= K, column j
*>          of the matrix A was interchanged with column JPIV(j).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors.
*> \endverbatim
*>
*> \param[in,out] VN1
*> \verbatim
*>          VN1 is DOUBLE PRECISION array, dimension (N)
*>          The vector with the partial column norms.
*> \endverbatim
*>
*> \param[in,out] VN2
*> \verbatim
*>          VN2 is DOUBLE PRECISION array, dimension (N)
*>          The vector with the exact column norms.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*>          Used in DLARF subroutine to apply elementary
*>          reflector.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup laqp2rk
*
*> \par References:
*  ================
*> [1] A Level 3 BLAS QR factorization algorithm with column pivoting developed in 1996.
*> G. Quintana-Orti, Depto. de Informatica, Universidad Jaime I, Spain.
*> X. Sun, Computer Science Dept., Duke University, USA.
*> C. H. Bischof, Math. and Comp. Sci. Div., Argonne National Lab, USA.
*> A BLAS-3 version of the QR factorization with column pivoting.
*> LAPACK Working Note 114
*> \htmlonly
*> <a href="https://www.netlib.org/lapack/lawnspdf/lawn114.pdf">https://www.netlib.org/lapack/lawnspdf/lawn114.pdf</a>
*> \endhtmlonly
*> and in
*> SIAM J. Sci. Comput., 19(5):1486-1494, Sept. 1998.
*> \htmlonly
*> <a href="https://doi.org/10.1137/S1064827595296732">https://doi.org/10.1137/S1064827595296732</a>
*> \endhtmlonly
*>
*> [2] A partial column norm updating strategy developed in 2006.
*> Z. Drmac and Z. Bujanovic, Dept. of Math., University of Zagreb, Croatia.
*> On the failure of rank revealing QR factorization software – a case study.
*> LAPACK Working Note 176.
*> \htmlonly
*> <a href="http://www.netlib.org/lapack/lawnspdf/lawn176.pdf">http://www.netlib.org/lapack/lawnspdf/lawn176.pdf</a>
*> \endhtmlonly
*
*  =====================================================================
      SUBROUTINE DLAQP2RK( M, N, NRHS, IOFFSET, KMAX, ABSTOL, RELTOL,
     $                     KP1, MAXC2NRM, A, LDA, KF, MAXC2NRMK,
     $                     RELMAXC2NRMK, JPIV, TAU, VN1, VN2, WORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            IOFFSET, KP1, KF, KMAX, LDA, M, N, NRHS
      DOUBLE PRECISION   ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK,
     $                   RELTOL
*     ..
*     .. Array Arguments ..
      INTEGER            JPIV( * )
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), VN1( * ), VN2( * ),
     $                   WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ITEMP, J, JMAXC2NRM, K, KP, MINMNFACT,
     $                   MINMNUPDT
      DOUBLE PRECISION   AIK, TEMP, TEMP2, TOL3Z
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLARFG, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH, DNRM2
      EXTERNAL           IDAMAX, DLAMCH, DNRM2
*     ..
*     .. Executable Statements ..
*
*     MINMNFACT in the smallest dimension of the submatrix
*     A(IOFFSET+1:M,1:N) to be factorized.
*
      MINMNFACT = MIN( M-IOFFSET, N )
      MINMNUPDT = MIN( M-IOFFSET, N+NRHS )
      KMAX = MIN( KMAX, MINMNFACT )
      TOL3Z = SQRT( DLAMCH( 'Epsilon' ) )
*
*     Compute factorization.
*
      DO K = 1, KMAX
*
         I = IOFFSET + K
*
         IF( IOFFSET.EQ.0 .AND. K.EQ.1 ) THEN
*
*           If we are at the first column of the original whole matrix A.
*
            KP = KP1
            MAXC2NRMK = MAXC2NRM
            RELMAXC2NRMK = ONE
*
         ELSE
*
*           Determine the pivot column at K-th step, i.e. the index
*           of the column with the maximum 2-norm in the
*           submatrix A(I:M,K:N).
*
            KP = ( K-1 ) + IDAMAX( N-K+1, VN1( K ), 1 )
*
*           Determine the maximum column 2-norm and the relative maximum
*           column 2-norm of the submatrix A(I:M,K:N) at step K.
*
            MAXC2NRMK = VN1( KP )
            RELMAXC2NRMK =  MAXC2NRMK / MAXC2NRM
*
         END IF
*
*     ==================================================================
*
*        Test for the second and third stopping criteria.
*        NOTE: There is no need to test for ABSTOL.GE.ZERO, since
*        MAXC2NRMK is non-negative. Similarly, there is no need
*        to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is
*        non-negative.
*

         IF( MAXC2NRMK.LE.ABSTOL .OR. RELMAXC2NRMK.LE.RELTOL ) THEN
*
*           Exit the loop
*
            EXIT
         END IF
*
*     ==================================================================
*
*        If the pivot column is not the first column of the
*        subblock A(1:M,K:N):
*        1) swap the K-th column and the KP-th pivot column
*           in A(1:M,1:N);
*        2) copy the K-th element into the KP-th element of the partial
*           and exact 2-norm vectors VN1 and VN2. ( Swap is not needed
*           for VN1 and VN2 since we use the element with the index
*           larger than K in the next loop step.)
*        3) Save the pivot interchange with the indices relative to the
*           the original matrix A, not the block A(1:M,1:N).
*
         IF( KP.NE.K ) THEN
            CALL DSWAP( M, A( 1, KP ), 1, A( 1, K ), 1 )
            VN1( KP ) = VN1( K )
            VN2( KP ) = VN2( K )
            ITEMP = JPIV( KP )
            JPIV( KP ) = JPIV( K )
            JPIV( K ) = ITEMP
         END IF
*
*        Generate elementary reflector H(K) using the column A(I:M,K),
*        if the column has more than one element, otherwise
*        the elementary reflector would be an identity matrix,
*        and TAU(K) = ZERO.
*
         IF( K.LT.M ) THEN
            CALL DLARFG( M-I+1, A( I, K ), A( I+1, K ), 1,
     $                   TAU( K ) )
         ELSE
            TAU( K ) = ZERO
         END IF
*
*       Apply H(K)**T to A(I:M,K+1:N+NRHS) from the left.
*       ( If M >= N, then at K = N there is no residual matrix,
*         i.e. no columns of A to update, only columns of B )
*         If M < N, then at K = M-IOFFSET, I = M and we have a
*         one-row residual matrix in A and the elementary
*         reflector is a unit matrix, TAU(K) = ZERO, i.e. no update
*         is needed for the residual matrix in A and the
*         right-hand-side-matrix in B.
*         Therefore, we update only if
*         K < MINMNUPDT = min(M-IOFFSET, N+NRHS)
*         condition is satisfied, not only K < N+NRHS )
*
         IF( K.LT.MINMNUPDT ) THEN
            AIK = A( I, K )
            A( I, K ) = ONE
            CALL DLARF( 'Left', M-I+1, N+NRHS-K, A( I, K ), 1,
     $                  TAU( K ), A( I, K+1 ), LDA, WORK( 1 ) )
            A( I, K ) = AIK
         END IF
*
         IF( K.LT.MINMNFACT ) THEN
*
*           Update the partial column 2-norms for the residual matrix,
*           only if the residual matrix A(I+1:M,K+1:N) exists, i.e.
*           when K < min(M-IOFFSET, N).
*
            DO J = K + 1, N
               IF( VN1( J ).NE.ZERO ) THEN
*
*                 NOTE: The following lines follow from the analysis in
*                 Lapack Working Note 176.
*
                  TEMP = ONE - ( ABS( A( I, J ) ) / VN1( J ) )**2
                  TEMP = MAX( TEMP, ZERO )
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                  IF( TEMP2 .LE. TOL3Z ) THEN
*
*                    Compute the column 2-norm for the partial
*                    column A(I+1:M,J) by explicitly computing it,
*                    and store it in both partial 2-norm vector VN1
*                    and exact column 2-norm vector VN2.
*
                     VN1( J ) = DNRM2( M-I, A( I+1, J ), 1 )
                     VN2( J ) = VN1( J )
*
                  ELSE
*
*                    Update the column 2-norm for the partial
*                    column A(I+1:M,J) by removing one
*                    element A(I,J) and store it in partial
*                    2-norm vector VN1.
*
                     VN1( J ) = VN1( J )*SQRT( TEMP )
*
                  END IF
               END IF
            END DO
*
         END IF
*
*     End factorization loop
*
      END DO
*
*     Set the number of factorized columns
*
      KF = K - 1
*
      IF( KF.EQ.KMAX ) THEN
*
*        All KMAX columns were factorized, no ABSTOL or RELTOL triggered.
*

         IF( KF.LT.MINMNFACT ) THEN
            JMAXC2NRM = KF + IDAMAX( N-KF, VN1( KF+1 ), 1 )
            MAXC2NRMK = VN1( JMAXC2NRM )
*
            IF( KF.EQ.0 ) THEN
               RELMAXC2NRMK = ONE
            ELSE
               RELMAXC2NRMK = MAXC2NRMK / MAXC2NRM
            END IF
*
         ELSE
            MAXC2NRMK = ZERO
            RELMAXC2NRMK = ZERO
         END IF
*

*
      END IF
*
*     Set TAU(KF+1:MINMN) to ZERO.
*
      DO J = KF + 1, MINMNFACT
         TAU( J ) = ZERO
      END DO
*
      RETURN
*
*     End of DLAQP2RK
*
      END
