*> \brief \b DLAQP3RK computes a step of truncated QR factorization with column pivoting of a real m-by-n matrix A using Level 3 BLAS and overwrites m-by-nrhs matrix B with Q**T * B.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLAQP3RK + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqp3rk.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqp3rk.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqp3rk.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*      SUBROUTINE DLAQP3RK( M, N, NRHS, IOFFSET, NB, KMAX, ABSTOL,
*     $                     RELTOL, KP1, MAXC2NRM, A, LDA, KB, DONE,
*     $                     KF, MAXC2NRMK, RELMAXC2NRMK,
*     $                     JPIV, TAU, VN1, VN2, AUXV, F, LDF, IWORK )
*      IMPLICIT NONE
*      LOGICAL            DONE
*      INTEGER            IOFFSET, KB, KF, KP1, LDA, LDF, M, KMAX, N,
*     $                   NB, NRHS
*      DOUBLE PRECISION   ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK,
*     $                   RELTOL
*
*     .. Scalar Arguments ..
*      LOGICAL            DONE
*      INTEGER            KB, LDA, LDF, M, N, NB, NRHS, IOFFSET
*      DOUBLE PRECISION   ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK,
*     $                   RELTOL
*     ..
*     .. Array Arguments ..
*      INTEGER            IWORK( * ), JPIV( * )
*      DOUBLE PRECISION   A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),
*     $                   VN1( * ), VN2( * )
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAQP3RK computes a step of truncated QR factorization with column
*> pivoting of a real M-by-N matrix A by using Level 3 BLAS.  The routine
*> tries to factorize NB columns from A starting from the row IOFFSET+1,
*> and updates all of the matrix with BLAS 3 xGEMM, the number of accually
*> factorized columns is returned in KB, KB <= NB.
*>
*> Cases when the number of factorized columns KB < NB:
*>
*> (1) In some cases, due to catastrophic cancellations, it cannot
*> factorize NB columns.  Hence, the actual number of factorized
*> columns is returned in KB.
*>
*> (2) Whenever the stopping criterion ABSTOL or RELTOL is satisfied,
*> the factorization is stopped, the logical DONE is returned
*> as TRUE.  The number of factorized columns which is smaller than NB
*> returned in KB.
*>
*> Block A(1:IOFFSET,1:N) is accordingly pivoted, but not factorized.
*>
*> The routine also overwrites the right-hand-sides B block stored
*> in A(IOFFSET+1:M,1:N+1:N+NRHS) with Q(K)**T * B.
*>
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
*>          The number of columns of the matrix A. N >= 0
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
*>          but no factorized. IOFFSET >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The number of columns to factorize.
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
*>          1. The elements in block A(IOFFSET+1:M,1:KB) below
*>             the diagonal,together with the array TAU, represent
*>             the orthogonal matrix Q(K) as a product of elementary
*>             reflectors.
*>          2. The block of the matrix A stored in A(IOFFSET+1:M,1:KB)
*>             is the triangular factor obtained.
*>          3. The block of the the matrix A stored in A(1:IOFFSET,1:N)
*>             has been accordingly pivoted, but no factorized.
*>          4. The rest of the array A, block A(IOFFSET+1:M,KB+1:N+NRHS).
*>             The left part A(IOFFSET+1:M,KB+1:N) of
*>             this block contains the residual of the matrix A, and
*>             the right part of the block A(IOFFSET+1:M,N+1:N+NRHS)
*>             contains the block of the right-hand-side matrix B. Both
*>             these blocks have been updated by multiplication from
*>             the left by Q**T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out]
*>
*> \verbatim
*>          DONE is LOGICAL
*>          TRUE, if the factorization completed,
*>          FALSE, otherwise.
*> \endverbatim
*
*> \param[out] KB
*> \verbatim
*>          KB is INTEGER
*>          The number of columns actually factorized.
*> \endverbatim
*>
*> \param[out] KF
*> \verbatim
*>          KF is INTEGER
*>          The number of columns of the original whole matrix A
*>          factorized.
*> \endverbatim
*
*> \param[out] MAXC2NRMK
*> \verbatim
*>          MAXC2NRMK is DOUBLE PRECISION
*>          The maximum column 2-norm of the residual matrix A22(K),
*>          when factorization stopped at rank K. MAXC2NRMK >= 0.
*>          ( Rank K is with respect to the original matrix A )
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
*>          Column pivot indices, for 1 <= j <= N, column j
*>          of the matrix A was interchanged with column JPIV(j).
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (NB)
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
*> \param[out] AUXV
*> \verbatim
*>          AUXV is DOUBLE PRECISION array, dimension (NB)
*>          Auxiliary vector.
*> \endverbatim
*>
*> \param[out] F
*> \verbatim
*>          F is DOUBLE PRECISION array, dimension (LDF,NB)
*>          Matrix F**T = L*Y**T*A.
*> \endverbatim
*>
*> \param[in] LDF
*> \verbatim
*>          LDF is INTEGER
*>          The leading dimension of the array F. LDF >= max(1,N+NRHS).
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N-1).
*>          Is a work array. ( IWORK is used to store indices
*>          of "bad" columns for norm downdating in the residual
*>          matrix ).
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
*> \ingroup laqp3rk
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
      SUBROUTINE DLAQP3RK( M, N, NRHS, IOFFSET, NB, KMAX, ABSTOL,
     $                     RELTOL, KP1, MAXC2NRM, A, LDA, DONE, KB,
     $                     KF, MAXC2NRMK, RELMAXC2NRMK,
     $                     JPIV, TAU, VN1, VN2, AUXV, F, LDF, IWORK )
      IMPLICIT NONE
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            DONE
      INTEGER            IOFFSET, KB, KF, KP1, LDA, LDF, M, KMAX, N,
     $                   NB, NRHS
      DOUBLE PRECISION   ABSTOL, MAXC2NRM, MAXC2NRMK, RELMAXC2NRMK,
     $                   RELTOL
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), JPIV( * )
      DOUBLE PRECISION   A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ),
     $                   VN1( * ), VN2( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ITEMP, J, K, MINMNFACT, MINMNUPDT,
     $                   LSTICC, KP, I, IF
      DOUBLE PRECISION   AIK, TEMP, TEMP2, TOL3Z
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DLARFG, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
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
      TOL3Z = SQRT( DLAMCH( 'Epsilon' ) )
*
*     Compute factorization in a while loop over NB columns,
*     K is the column index in the block A(1:M,1:N).
*
      K = 0
      LSTICC = 0
      DONE = .FALSE.
*
      DO WHILE ( K.LT.NB .AND. LSTICC.EQ.0 )
         K = K + 1
         I = IOFFSET + K
*
         IF( I.EQ.1 ) THEN
*
*           We are at the first column of the original whole matrix A,
*           therefore we use the computed KP1 and MAXC2NRM from the
*           main routine.
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
*           TODO: optimize RELMAXC2NRMK
            RELMAXC2NRMK =  MAXC2NRMK / MAXC2NRM
*
         END IF
*
*     ==================================================================
*
*        Quick return, if the submatrix A(I:M,K:N) is
*        a zero matrix. We need to check it only if the column index
*        (same as row index) is larger than 2, since the condition for
*        the whole original matrix is checked in the main routine.
*
         IF( I.NE.1 .AND. MAXC2NRMK.EQ.ZERO ) THEN


         WRITE(*,*) "$$$$$$ DLAQP3RK zero submatrix, IOFFSET, K= ",
     $                     IOFFSET, K
*
            DONE = .TRUE.
*
*           Set KB, the number of factorized columns in the block;
*           Set IF, the number of processed rows in the block, which is
*           the same as the number of rows in the original whole
*           matrix A;
*           Set KF, the number of factorized columns in the original
*           whole matrix A.
*           TODO: fix USETOL
            IF( MAXC2NRMK.LE.ABSTOL .OR. RELMAXC2NRMK.LE.RELTOL ) THEN


               WRITE(*,*)
     $         "$$$$$$$$ DLAQP3RK zero submatrix (ABSTOL, K)= ",
     $         ABSTOL,  K
*
               KB = K - 1
               IF = I - 1
               KF = IOFFSET + KB
*
            ELSE
*
               KB = K - 1
               IF = I - 1
               KF = KMAX
*
            END IF
*
*           There is no need to apply the block reflector to the
*           residual of the matrix A stored in A(KB+1:M,KB+1:N), since
*           the submatrix is zero and we stop the computation.  But,
*           we need to apply the block reflector to the residual right
*           hand sides stored in A(KB+1:M,N+1:N+NRHS), if the residual
*           right hand sides exist.  This happens
*           when ( NRHS != 0 AND KB <= (M-IOFFSET) ):
*
*           A(I+1:M,N+1:N+NRHS) := A(I+1:M,N+1:N+NRHS) -
*                            A(I+1:M,1:KB) * F(N+1:N+NRHS,1:KB)**T.
*
            IF( NRHS.GT.0 .AND. KB.LT.(M-IOFFSET) ) THEN


            WRITE(*,*) "$$$$$$$$$$ DLAQP3RK block reflector ",
     $              "(M-IF, NRHS, KB)", M-IF, NRHS, KB

*
               CALL DGEMM( 'No transpose', 'Transpose', M-IF, NRHS,
     $                      KB, -ONE, A( IF+1, 1 ), LDA, F( N+1, 1 ),
     $                      LDF, ONE, A( IF+1, N+1 ), LDA )
            END IF
*
*           There is no need to recompute the 2-norm of the
*           difficult columns, since we stop the factorization.
*
*           Set TAUs corresponding to the columns that were not
*           factorized to ZERO, i.e. set TAU(KB+1:MINMNFACT) to ZERO,
*           which is equivalent to seting TAU(K:MINMNFACT) to ZERO.
*
            DO J = K, MINMNFACT
               TAU( J ) = ZERO
            END DO
*
*           Return from the routine.
*
            RETURN
*
         END IF
*
*     ==================================================================
*
*        Test for the second and third tolerance stopping criteria.
*        NOTE: There is no need to test for ABSTOL.GE.ZERO, since
*        MAXC2NRMK is non-negative. Similarly, there is no need
*        to test for RELTOL.GE.ZERO, since RELMAXC2NRMK is
*        non-negative.
*
         IF( MAXC2NRMK.LE.ABSTOL .OR. RELMAXC2NRMK.LE.RELTOL ) THEN
*
            DONE = .TRUE.
*
*           Set the number of factorized columns in the block.
*
            K = K - 1
*
*           Exit the loop.
*           After the loop, there is a code:
*              1) to apply the block reflector via GEMM to the residual
*                 of the matrix A and to the right hand sides B;
*              2) to recompute the 2-norm of the difficult columns;
*              3) to zero out the remaining TAUs.
*           TODO: change exit??
            EXIT
*
         END IF
*
*     ==================================================================
*
*        If the pivot column is not the first column of the
*        subblock A(1:M,K:N):
*        1) swap the K-th column and the KP-th pivot column
*           in A(1:M,1:N);
*        2) swap the K-th row and the KP-th row in F(1:N,1:K-1)
*        3) copy the K-th element into the KP-th element of the partial
*           and exact 2-norm vectors VN1 and VN2. (Swap is not needed
*           for VN1 and VN2 since we use the element with the index
*           larger than K in the next loop step.)
*        4) Save the pivot interchange with the indices relative to the
*           the original matrix A, not the block A(1:M,1:N).
*
         IF( KP.NE.K ) THEN
            CALL DSWAP( M, A( 1, KP ), 1, A( 1, K ), 1 )
            CALL DSWAP( K-1, F( KP, 1 ), LDF, F( K, 1 ), LDF )
            VN1( KP ) = VN1( K )
            VN2( KP ) = VN2( K )
            ITEMP = JPIV( KP )
            JPIV( KP ) = JPIV( K )
            JPIV( K ) = ITEMP
         END IF
*
*        Apply previous Householder reflectors to column K:
*        A(I:M,K) := A(I:M,K) - A(I:M,1:K-1)*F(K,1:K-1)**T.
*
         IF( K.GT.1 ) THEN
            CALL DGEMV( 'No transpose', M-I+1, K-1, -ONE, A( I, 1 ),
     $                  LDA, F( K, 1 ), LDF, ONE, A( I, K ), 1 )
         END IF
*
*        Generate elementary reflector H(k) using the column A(I:M,K).
*
         IF( I.LT.M ) THEN
            CALL DLARFG( M-I+1, A( I, K ), A( I+1, K ), 1, TAU( K ) )
         ELSE
            TAU( K ) = ZERO
         END IF
*
         AIK = A( I, K )
         A( I, K ) = ONE
*        ===============================================================
*
*        Compute the current K-th column of F:
*          1) F(K+1:N,K) := tau(K) * A(I:M,K+1:N)**T * A(I:M,K).
*
         IF( K.LT.N+NRHS ) THEN
            CALL DGEMV( 'Transpose', M-I+1, N+NRHS-K, TAU( K ),
     $                  A( I, K+1 ), LDA, A( I, K ), 1, ZERO,
     $                  F( K+1, K ), 1 )
         END IF
*
*           2) Zero out elements above and on the diagonal of the
*              column K in matrix F, i.e elements F(1:K,K).
*
         DO J = 1, K
            F( J, K ) = ZERO
         END DO
*
*         3) Incremental updating of the K-th column of F:
*        F(1:N,K) := F(1:N,K) - tau(K) * F(1:N,1:K-1) * A(I:M,1:K-1)**T
*                    * A(I:M,K).
*
         IF( K.GT.1 ) THEN
            CALL DGEMV( 'Transpose', M-I+1, K-1, -TAU( K ), A( I, 1 ),
     $                  LDA, A( I, K ), 1, ZERO, AUXV( 1 ), 1 )
*
            CALL DGEMV( 'No transpose', N+NRHS, K-1, ONE,
     $                  F( 1, 1 ), LDF, AUXV( 1 ), 1, ONE,
     $                  F( 1, K ), 1 )
         END IF
*
*        ===============================================================
*
*        Update the current I-th row of A:
*        A(I,K+1:N) := A(I,K+1:N) - A(I,1:K)*F(K+1:N,1:K)**T.
*
         IF( K.LT.N+NRHS ) THEN
            CALL DGEMV( 'No transpose', N+NRHS-K, K, -ONE,
     $                  F( K+1, 1 ), LDF, A( I, 1 ), LDA, ONE,
     $                  A( I, K+1 ), LDA )
         END IF
*
         A( I, K ) = AIK
*
*        Update the partial column 2-norms for the residual matrix,
*        only if the residual matrix A(I+1:M,K+1:N) exists, i.e.
*        when K < MINMNFACT = min( M-IOFFSET, N ).
*
         IF( K.LT.MINMNFACT ) THEN
*
            DO J = K + 1, N
               IF( VN1( J ).NE.ZERO ) THEN
*
*                 NOTE: The following lines follow from the analysis in
*                 Lapack Working Note 176.
*
                  TEMP = ABS( A( I, J ) ) / VN1( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                  IF( TEMP2.LE.TOL3Z ) THEN
*
*                    At J-index, we have a difficult column for the
*                    update of the 2-norm. Save the index of the previous
*                    difficult column in IWORK(J-1).
*                    NOTE: ILSTCC > 1, threfore we can use IWORK only
*                    with N-1 elements, where the elements are
*                    shifted by 1 to the left.
*
                     IWORK( J-1 ) = LSTICC
*
*                    Set the index of the last difficult column LSTICC.
*
                     LSTICC = J
*
                  ELSE
                     VN1( J ) = VN1( J )*SQRT( TEMP )
                  END IF
               END IF
            END DO
*
         END IF
*
*        End of while loop.
*
      END DO
*
*     Now, afler the loop:
*        KB is the number of factorized columns in the block,
*        KF is the number or factorized columns in the original
*           whole matrix A,
*        I  is the number of processed rows in the block which is
*           the same as the the numerb of processed rows in
*           the original whole matrix A.
*
      KB = K
      KF = IOFFSET + KB
      I = IOFFSET + KB

*     Apply the block reflector to the residual of the matrix A
*     and right hand sides B, if the residual matrix and
*     and/or the residual right hand sides exist, i.e.
*     if the submatrix A(I+1:M,KB+1:N+NRHS) exists. This happens
*     when KB < MINMNUPDT = min( M-IOFFSET, N+NRHS ):
*
*     A(I+1:M,K+1:N+NRHS) := A(I+1:M,KB+1:N+NRHS) -
*                         A(I+1:M,1:KB) * F(KB+1:N+NRHS,1:KB)**T.
*
      IF( KB.LT.MINMNUPDT ) THEN
*
         CALL DGEMM( 'No transpose', 'Transpose', M-I, N+NRHS-KB, KB,
     $               -ONE, A( I+1, 1 ), LDA, F( KB+1, 1 ), LDF, ONE,
     $               A( I+1, KB+1 ), LDA )
      END IF
*
*     Recompute the 2-norm of the difficult columns.
*     Loop over the index of the difficult columns from the largest
*     to the smallest index.
*
      DO WHILE( LSTICC.GT.0 )
*
*        LSTICC is the index of the last difficult column is greater
*        than 1.
*        ITEMP is the index of the previous difficult column.
*
         ITEMP = IWORK( LSTICC-1 )
*
*        Compute the 2-norm explicilty for the last difficult column and
*        save it in the partial and exact 2-norm vectors VN1 and VN2.
*
*        NOTE: The computation of VN1( LSTICC ) relies on the fact that
*        DNRM2 does not fail on vectors with norm below the value of
*        SQRT(DLAMCH('S'))
*
         VN1( LSTICC ) = DNRM2( M-I, A( I+1, LSTICC ), 1 )
         VN2( LSTICC ) = VN1( LSTICC )
*
*        Downdate the index of the last difficult column to
*        the index of the previous difficult column.
*
         LSTICC = ITEMP
      END DO
*
*     If done, set TAU(KB+1:MINMNFACT) to ZERO.
*
      IF( DONE ) THEN
         DO J = KB + 1, MINMNFACT
            TAU( J ) = ZERO
         END DO
      END IF
*
      RETURN
*
*     End of DLAQP3RK
*
      END
