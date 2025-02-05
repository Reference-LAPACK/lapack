*> \brief <b> DGGLSE solves overdetermined or underdetermined systems for OTHER matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DGGLSE + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgglse.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgglse.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgglse.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK,
*                          INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( * ), D( * ),
*      $                   WORK( * ), X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGGLSE solves the linear equality-constrained least squares (LSE)
*> problem:
*>
*>         minimize || c - A*x ||_2   subject to   B*x = d
*>
*> where A is an M-by-N matrix, B is a P-by-N matrix, c is a given
*> M-vector, and d is a given P-vector. It is assumed that
*> P <= N <= M+P, and
*>
*>          rank(B) = P and  rank( (A) ) = N.
*>                               ( (B) )
*>
*> These conditions ensure that the LSE problem has a unique solution,
*> which is obtained using a generalized RQ factorization of the
*> matrices (B, A) given by
*>
*>    B = (0 R)*Q,   A = Z*T*Q.
*>
*> Callers of this subroutine should note that the singularity/rank-deficiency checks
*> implemented in this subroutine are rudimentary. The DTRTRS subroutine called by this
*> subroutine only signals a failure due to singularity if the problem is exactly singular.
*>
*> It is conceivable for one (or more) of the factors involved in the generalized RQ
*> factorization of the pair (B, A) to be subnormally close to singularity without this
*> subroutine signalling an error. The solutions computed for such almost-rank-deficient
*> problems may be less accurate due to a loss of numerical precision.
*> 
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrices A and B. N >= 0.
*> \endverbatim
*>
*> \param[in] P
*> \verbatim
*>          P is INTEGER
*>          The number of rows of the matrix B. 0 <= P <= N <= M+P.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the M-by-N matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(M,N)-by-N upper trapezoidal matrix T.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,N)
*>          On entry, the P-by-N matrix B.
*>          On exit, the upper triangle of the subarray B(1:P,N-P+1:N)
*>          contains the P-by-P upper triangular matrix R.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B. LDB >= max(1,P).
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (M)
*>          On entry, C contains the right hand side vector for the
*>          least squares part of the LSE problem.
*>          On exit, the residual sum of squares for the solution
*>          is given by the sum of squares of elements N-P+1 to M of
*>          vector C.
*> \endverbatim
*>
*> \param[in,out] D
*> \verbatim
*>          D is DOUBLE PRECISION array, dimension (P)
*>          On entry, D contains the right hand side vector for the
*>          constrained equation.
*>          On exit, D is destroyed.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (N)
*>          On exit, X is the solution of the LSE problem.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,M+N+P).
*>          For optimum performance LWORK >= P+min(M,N)+max(M,N)*NB,
*>          where NB is an upper bound for the optimal blocksizes for
*>          DGEQRF, SGERQF, DORMQR and SORMRQ.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*>          = 1:  the upper triangular factor R associated with B in the
*>                generalized RQ factorization of the pair (B, A) is exactly
*>                singular, so that rank(B) < P; the least squares
*>                solution could not be computed.
*>          = 2:  the (N-P) by (N-P) part of the upper trapezoidal factor
*>                T associated with A in the generalized RQ factorization
*>                of the pair (B, A) is exactly singular, so that
*>                rank( (A) ) < N; the least squares solution could not
*>                    ( (B) )
*>                be computed.
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
*> \ingroup gglse
*
*  =====================================================================
      SUBROUTINE DGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK,
     $                   LWORK,
     $                   INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LDB, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( * ), D( * ),
     $                   WORK( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3,
     $                   NB4, NR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DGEMV, DGGRQF, DORMQR,
     $                   DORMRQ,
     $                   DTRMV, DTRTRS, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          INT, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( P.LT.0 .OR. P.GT.N .OR. P.LT.N-M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -7
      END IF
*
*     Calculate workspace
*
      IF( INFO.EQ.0) THEN
         IF( N.EQ.0 ) THEN
            LWKMIN = 1
            LWKOPT = 1
         ELSE
            NB1 = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 )
            NB3 = ILAENV( 1, 'DORMQR', ' ', M, N, P, -1 )
            NB4 = ILAENV( 1, 'DORMRQ', ' ', M, N, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = P + MN + MAX( M, N )*NB
         END IF
         WORK( 1 ) = LWKOPT
*
         IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGLSE', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Compute the GRQ factorization of matrices B and A:
*
*            B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P
*                        N-P  P                     (  0  R22 ) M+P-N
*                                                      N-P  P
*
*     where T12 and R11 are upper triangular, and Q and Z are
*     orthogonal.
*
      CALL DGGRQF( P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ),
     $             WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = INT( WORK( P+MN+1 ) )
*
*     Update c = Z**T *c = ( c1 ) N-P
*                          ( c2 ) M+P-N
*
      CALL DORMQR( 'Left', 'Transpose', M, 1, MN, A, LDA,
     $             WORK( P+1 ),
     $             C, MAX( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO )
      LOPT = MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
*     Solve T12*x2 = d for x2
*
      IF( P.GT.0 ) THEN
         CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', P, 1,
     $                B( 1, N-P+1 ), LDB, D, P, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 1
            RETURN
         END IF
*
*        Put the solution in X
*
         CALL DCOPY( P, D, 1, X( N-P+1 ), 1 )
*
*        Update c1
*
         CALL DGEMV( 'No transpose', N-P, P, -ONE, A( 1, N-P+1 ),
     $               LDA,
     $               D, 1, ONE, C, 1 )
      END IF
*
*     Solve R11*x1 = c1 for x1
*
      IF( N.GT.P ) THEN
         CALL DTRTRS( 'Upper', 'No transpose', 'Non-unit', N-P, 1,
     $                A, LDA, C, N-P, INFO )
*
         IF( INFO.GT.0 ) THEN
            INFO = 2
            RETURN
         END IF
*
*        Put the solutions in X
*
         CALL DCOPY( N-P, C, 1, X, 1 )
      END IF
*
*     Compute the residual vector:
*
      IF( M.LT.N ) THEN
         NR = M + P - N
         IF( NR.GT.0 )
     $      CALL DGEMV( 'No transpose', NR, N-M, -ONE, A( N-P+1,
     $                  M+1 ),
     $                  LDA, D( NR+1 ), 1, ONE, C( N-P+1 ), 1 )
      ELSE
         NR = P
      END IF
      IF( NR.GT.0 ) THEN
         CALL DTRMV( 'Upper', 'No transpose', 'Non unit', NR,
     $               A( N-P+1, N-P+1 ), LDA, D, 1 )
         CALL DAXPY( NR, -ONE, D, 1, C( N-P+1 ), 1 )
      END IF
*
*     Backward transformation x = Q**T*x
*
      CALL DORMRQ( 'Left', 'Transpose', N, 1, P, B, LDB, WORK( 1 ),
     $             X,
     $             N, WORK( P+MN+1 ), LWORK-P-MN, INFO )
      WORK( 1 ) = P + MN + MAX( LOPT, INT( WORK( P+MN+1 ) ) )
*
      RETURN
*
*     End of DGGLSE
*
      END
