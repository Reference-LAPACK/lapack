      SUBROUTINE DGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV,
     $                    R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS,
     $                    ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,
     $                    WORK, IWORK, INFO )
*
*     -- LAPACK routine (version 3.2.1)                                 --
*     -- Contributed by James Demmel, Deaglan Halligan, Yozo Hida and --
*     -- Jason Riedy of Univ. of California Berkeley.                 --
*     -- April 2009                                                   --
*
*     -- LAPACK is a software package provided by Univ. of Tennessee, --
*     -- Univ. of California Berkeley and NAG Ltd.                    --
*
      IMPLICIT NONE
*     ..
*     .. Scalar Arguments ..
      CHARACTER          TRANS, EQUED
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,
     $                   N_ERR_BNDS
      DOUBLE PRECISION   RCOND
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   X( LDX , * ), WORK( * )
      DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),
     $                   ERR_BNDS_NORM( NRHS, * ),
     $                   ERR_BNDS_COMP( NRHS, * )
*     ..
*
*     Purpose
*     =======
*
*     DGERFSX improves the computed solution to a system of linear
*     equations and provides error bounds and backward error estimates
*     for the solution.  In addition to normwise error bound, the code
*     provides maximum componentwise error bound if possible.  See
*     comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the
*     error bounds.
*
*     The original system of linear equations may have been equilibrated
*     before calling this routine, as described by arguments EQUED, R
*     and C below. In this case, the solution and error bounds returned
*     are for the original unequilibrated system.
*
*     Arguments
*     =========
*
*     Some optional parameters are bundled in the PARAMS array.  These
*     settings determine how refinement is performed, but often the
*     defaults are acceptable.  If the defaults are acceptable, users
*     can pass NPARAMS = 0 which prevents the source code from accessing
*     the PARAMS argument.
*
*     TRANS   (input) CHARACTER*1
*     Specifies the form of the system of equations:
*       = 'N':  A * X = B     (No transpose)
*       = 'T':  A**T * X = B  (Transpose)
*       = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*     EQUED   (input) CHARACTER*1
*     Specifies the form of equilibration that was done to A
*     before calling this routine. This is needed to compute
*     the solution and error bounds correctly.
*       = 'N':  No equilibration
*       = 'R':  Row equilibration, i.e., A has been premultiplied by
*               diag(R).
*       = 'C':  Column equilibration, i.e., A has been postmultiplied
*               by diag(C).
*       = 'B':  Both row and column equilibration, i.e., A has been
*               replaced by diag(R) * A * diag(C).
*               The right hand side B has been changed accordingly.
*
*     N       (input) INTEGER
*     The order of the matrix A.  N >= 0.
*
*     NRHS    (input) INTEGER
*     The number of right hand sides, i.e., the number of columns
*     of the matrices B and X.  NRHS >= 0.
*
*     A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*     The original N-by-N matrix A.
*
*     LDA     (input) INTEGER
*     The leading dimension of the array A.  LDA >= max(1,N).
*
*     AF      (input) DOUBLE PRECISION array, dimension (LDAF,N)
*     The factors L and U from the factorization A = P*L*U
*     as computed by DGETRF.
*
*     LDAF    (input) INTEGER
*     The leading dimension of the array AF.  LDAF >= max(1,N).
*
*     IPIV    (input) INTEGER array, dimension (N)
*     The pivot indices from DGETRF; for 1<=i<=N, row i of the
*     matrix was interchanged with row IPIV(i).
*
*     R       (input or output) DOUBLE PRECISION array, dimension (N)
*     The row scale factors for A.  If EQUED = 'R' or 'B', A is
*     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
*     is not accessed.  R is an input argument if FACT = 'F';
*     otherwise, R is an output argument.  If FACT = 'F' and
*     EQUED = 'R' or 'B', each element of R must be positive.
*     If R is output, each element of R is a power of the radix.
*     If R is input, each element of R should be a power of the radix
*     to ensure a reliable solution and error estimates. Scaling by
*     powers of the radix does not cause rounding errors unless the
*     result underflows or overflows. Rounding errors during scaling
*     lead to refining with a matrix that is not equivalent to the
*     input matrix, producing error estimates that may not be
*     reliable.
*
*     C       (input or output) DOUBLE PRECISION array, dimension (N)
*     The column scale factors for A.  If EQUED = 'C' or 'B', A is
*     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
*     is not accessed.  C is an input argument if FACT = 'F';
*     otherwise, C is an output argument.  If FACT = 'F' and
*     EQUED = 'C' or 'B', each element of C must be positive.
*     If C is output, each element of C is a power of the radix.
*     If C is input, each element of C should be a power of the radix
*     to ensure a reliable solution and error estimates. Scaling by
*     powers of the radix does not cause rounding errors unless the
*     result underflows or overflows. Rounding errors during scaling
*     lead to refining with a matrix that is not equivalent to the
*     input matrix, producing error estimates that may not be
*     reliable.
*
*     B       (input) DOUBLE PRECISION array, dimension (LDB,NRHS)
*     The right hand side matrix B.
*
*     LDB     (input) INTEGER
*     The leading dimension of the array B.  LDB >= max(1,N).
*
*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,NRHS)
*     On entry, the solution matrix X, as computed by DGETRS.
*     On exit, the improved solution matrix X.
*
*     LDX     (input) INTEGER
*     The leading dimension of the array X.  LDX >= max(1,N).
*
*     RCOND   (output) DOUBLE PRECISION
*     Reciprocal scaled condition number.  This is an estimate of the
*     reciprocal Skeel condition number of the matrix A after
*     equilibration (if done).  If this is less than the machine
*     precision (in particular, if it is zero), the matrix is singular
*     to working precision.  Note that the error may still be small even
*     if this number is very small and the matrix appears ill-
*     conditioned.
*
*     BERR    (output) DOUBLE PRECISION array, dimension (NRHS)
*     Componentwise relative backward error.  This is the
*     componentwise relative backward error of each solution vector X(j)
*     (i.e., the smallest relative change in any element of A or B that
*     makes X(j) an exact solution).
*
*     N_ERR_BNDS (input) INTEGER
*     Number of error bounds to return for each right hand side
*     and each type (normwise or componentwise).  See ERR_BNDS_NORM and
*     ERR_BNDS_COMP below.
*
*     ERR_BNDS_NORM  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)
*     For each right-hand side, this array contains information about
*     various error bounds and condition numbers corresponding to the
*     normwise relative error, which is defined as follows:
*
*     Normwise relative error in the ith solution vector:
*             max_j (abs(XTRUE(j,i) - X(j,i)))
*            ------------------------------
*                  max_j abs(X(j,i))
*
*     The array is indexed by the type of error information as described
*     below. There currently are up to three pieces of information
*     returned.
*
*     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
*     right-hand side.
*
*     The second index in ERR_BNDS_NORM(:,err) contains the following
*     three fields:
*     err = 1 "Trust/don't trust" boolean. Trust the answer if the
*              reciprocal condition number is less than the threshold
*              sqrt(n) * dlamch('Epsilon').
*
*     err = 2 "Guaranteed" error bound: The estimated forward error,
*              almost certainly within a factor of 10 of the true error
*              so long as the next entry is greater than the threshold
*              sqrt(n) * dlamch('Epsilon'). This error bound should only
*              be trusted if the previous boolean is true.
*
*     err = 3  Reciprocal condition number: Estimated normwise
*              reciprocal condition number.  Compared with the threshold
*              sqrt(n) * dlamch('Epsilon') to determine if the error
*              estimate is "guaranteed". These reciprocal condition
*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
*              appropriately scaled matrix Z.
*              Let Z = S*A, where S scales each row by a power of the
*              radix so all absolute row sums of Z are approximately 1.
*
*     See Lapack Working Note 165 for further details and extra
*     cautions.
*
*     ERR_BNDS_COMP  (output) DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)
*     For each right-hand side, this array contains information about
*     various error bounds and condition numbers corresponding to the
*     componentwise relative error, which is defined as follows:
*
*     Componentwise relative error in the ith solution vector:
*                    abs(XTRUE(j,i) - X(j,i))
*             max_j ----------------------
*                         abs(X(j,i))
*
*     The array is indexed by the right-hand side i (on which the
*     componentwise relative error depends), and the type of error
*     information as described below. There currently are up to three
*     pieces of information returned for each right-hand side. If
*     componentwise accuracy is not requested (PARAMS(3) = 0.0), then
*     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS .LT. 3, then at most
*     the first (:,N_ERR_BNDS) entries are returned.
*
*     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
*     right-hand side.
*
*     The second index in ERR_BNDS_COMP(:,err) contains the following
*     three fields:
*     err = 1 "Trust/don't trust" boolean. Trust the answer if the
*              reciprocal condition number is less than the threshold
*              sqrt(n) * dlamch('Epsilon').
*
*     err = 2 "Guaranteed" error bound: The estimated forward error,
*              almost certainly within a factor of 10 of the true error
*              so long as the next entry is greater than the threshold
*              sqrt(n) * dlamch('Epsilon'). This error bound should only
*              be trusted if the previous boolean is true.
*
*     err = 3  Reciprocal condition number: Estimated componentwise
*              reciprocal condition number.  Compared with the threshold
*              sqrt(n) * dlamch('Epsilon') to determine if the error
*              estimate is "guaranteed". These reciprocal condition
*              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
*              appropriately scaled matrix Z.
*              Let Z = S*(A*diag(x)), where x is the solution for the
*              current right-hand side and S scales each row of
*              A*diag(x) by a power of the radix so all absolute row
*              sums of Z are approximately 1.
*
*     See Lapack Working Note 165 for further details and extra
*     cautions.
*
*     NPARAMS (input) INTEGER
*     Specifies the number of parameters set in PARAMS.  If .LE. 0, the
*     PARAMS array is never referenced and default values are used.
*
*     PARAMS  (input / output) DOUBLE PRECISION array, dimension NPARAMS
*     Specifies algorithm parameters.  If an entry is .LT. 0.0, then
*     that entry will be filled with default value used for that
*     parameter.  Only positions up to NPARAMS are accessed; defaults
*     are used for higher-numbered parameters.
*
*       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
*            refinement or not.
*         Default: 1.0D+0
*            = 0.0 : No refinement is performed, and no error bounds are
*                    computed.
*            = 1.0 : Use the double-precision refinement algorithm,
*                    possibly with doubled-single computations if the
*                    compilation environment does not support DOUBLE
*                    PRECISION.
*              (other values are reserved for future use)
*
*       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
*            computations allowed for refinement.
*         Default: 10
*         Aggressive: Set to 100 to permit convergence using approximate
*                     factorizations or factorizations other than LU. If
*                     the factorization uses a technique other than
*                     Gaussian elimination, the guarantees in
*                     err_bnds_norm and err_bnds_comp may no longer be
*                     trustworthy.
*
*       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
*            will attempt to find a solution with small componentwise
*            relative error in the double-precision algorithm.  Positive
*            is true, 0.0 is false.
*         Default: 1.0 (attempt componentwise convergence)
*
*     WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*     IWORK   (workspace) INTEGER array, dimension (N)
*
*     INFO    (output) INTEGER
*       = 0:  Successful exit. The solution to every right-hand side is
*         guaranteed.
*       < 0:  If INFO = -i, the i-th argument had an illegal value
*       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
*         has been completed, but the factor U is exactly singular, so
*         the solution and error bounds could not be computed. RCOND = 0
*         is returned.
*       = N+J: The solution corresponding to the Jth right-hand side is
*         not guaranteed. The solutions corresponding to other right-
*         hand sides K with K > J may not be guaranteed as well, but
*         only the first such right-hand side is reported. If a small
*         componentwise error is not requested (PARAMS(3) = 0.0) then
*         the Jth right-hand side is the first with a normwise error
*         bound that is not guaranteed (the smallest J such
*         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
*         the Jth right-hand side is the first with either a normwise or
*         componentwise error bound that is not guaranteed (the smallest
*         J such that either ERR_BNDS_NORM(J,1) = 0.0 or
*         ERR_BNDS_COMP(J,1) = 0.0). See the definition of
*         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
*         about all of the right-hand sides check ERR_BNDS_NORM or
*         ERR_BNDS_COMP.
*
*     ==================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ITREF_DEFAULT, ITHRESH_DEFAULT
      DOUBLE PRECISION   COMPONENTWISE_DEFAULT, RTHRESH_DEFAULT
      DOUBLE PRECISION   DZTHRESH_DEFAULT
      PARAMETER          ( ITREF_DEFAULT = 1.0D+0 )
      PARAMETER          ( ITHRESH_DEFAULT = 10.0D+0 )
      PARAMETER          ( COMPONENTWISE_DEFAULT = 1.0D+0 )
      PARAMETER          ( RTHRESH_DEFAULT = 0.5D+0 )
      PARAMETER          ( DZTHRESH_DEFAULT = 0.25D+0 )
      INTEGER            LA_LINRX_ITREF_I, LA_LINRX_ITHRESH_I,
     $                   LA_LINRX_CWISE_I
      PARAMETER          ( LA_LINRX_ITREF_I = 1,
     $                   LA_LINRX_ITHRESH_I = 2 )
      PARAMETER          ( LA_LINRX_CWISE_I = 3 )
      INTEGER            LA_LINRX_TRUST_I, LA_LINRX_ERR_I,
     $                   LA_LINRX_RCOND_I
      PARAMETER          ( LA_LINRX_TRUST_I = 1, LA_LINRX_ERR_I = 2 )
      PARAMETER          ( LA_LINRX_RCOND_I = 3 )
*     ..
*     .. Local Scalars ..
      CHARACTER(1)       NORM
      LOGICAL            ROWEQU, COLEQU, NOTRAN
      INTEGER            J, TRANS_TYPE, PREC_TYPE, REF_TYPE
      INTEGER            N_NORMS
      DOUBLE PRECISION   ANORM, RCOND_TMP
      DOUBLE PRECISION   ILLRCOND_THRESH, ERR_LBND, CWISE_WRONG
      LOGICAL            IGNORE_CWISE
      INTEGER            ITHRESH
      DOUBLE PRECISION   RTHRESH, UNSTABLE_THRESH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DGECON, DLA_GERFSX_EXTENDED
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME, BLAS_FPINFO_X, ILATRANS, ILAPREC
      EXTERNAL           DLAMCH, DLANGE, DLA_GERCOND
      DOUBLE PRECISION   DLAMCH, DLANGE, DLA_GERCOND
      LOGICAL            LSAME
      INTEGER            BLAS_FPINFO_X
      INTEGER            ILATRANS, ILAPREC
*     ..
*     .. Executable Statements ..
*
*     Check the input parameters.
*
      INFO = 0
      TRANS_TYPE = ILATRANS( TRANS )
      REF_TYPE = INT( ITREF_DEFAULT )
      IF ( NPARAMS .GE. LA_LINRX_ITREF_I ) THEN
         IF ( PARAMS( LA_LINRX_ITREF_I ) .LT. 0.0D+0 ) THEN
            PARAMS( LA_LINRX_ITREF_I ) = ITREF_DEFAULT
         ELSE
            REF_TYPE = PARAMS( LA_LINRX_ITREF_I )
         END IF
      END IF
*
*     Set default parameters.
*
      ILLRCOND_THRESH = DBLE( N ) * DLAMCH( 'Epsilon' )
      ITHRESH = INT( ITHRESH_DEFAULT )
      RTHRESH = RTHRESH_DEFAULT
      UNSTABLE_THRESH = DZTHRESH_DEFAULT
      IGNORE_CWISE = COMPONENTWISE_DEFAULT .EQ. 0.0D+0
*
      IF ( NPARAMS.GE.LA_LINRX_ITHRESH_I ) THEN
         IF ( PARAMS( LA_LINRX_ITHRESH_I ).LT.0.0D+0 ) THEN
            PARAMS( LA_LINRX_ITHRESH_I ) = ITHRESH
         ELSE
            ITHRESH = INT( PARAMS( LA_LINRX_ITHRESH_I ) )
         END IF
      END IF
      IF ( NPARAMS.GE.LA_LINRX_CWISE_I ) THEN
         IF ( PARAMS( LA_LINRX_CWISE_I ).LT.0.0D+0 ) THEN
            IF ( IGNORE_CWISE ) THEN
               PARAMS( LA_LINRX_CWISE_I ) = 0.0D+0
            ELSE
               PARAMS( LA_LINRX_CWISE_I ) = 1.0D+0
            END IF
         ELSE
            IGNORE_CWISE = PARAMS( LA_LINRX_CWISE_I ) .EQ. 0.0D+0
         END IF
      END IF
      IF ( REF_TYPE .EQ. 0 .OR. N_ERR_BNDS .EQ. 0 ) THEN
         N_NORMS = 0
      ELSE IF ( IGNORE_CWISE ) THEN
         N_NORMS = 1
      ELSE
         N_NORMS = 2
      END IF
*
      NOTRAN = LSAME( TRANS, 'N' )
      ROWEQU = LSAME( EQUED, 'R' ) .OR. LSAME( EQUED, 'B' )
      COLEQU = LSAME( EQUED, 'C' ) .OR. LSAME( EQUED, 'B' )
*
*     Test input parameters.
*
      IF( TRANS_TYPE.EQ.-1 ) THEN
        INFO = -1
      ELSE IF( .NOT.ROWEQU .AND. .NOT.COLEQU .AND.
     $         .NOT.LSAME( EQUED, 'N' ) ) THEN
        INFO = -2
      ELSE IF( N.LT.0 ) THEN
        INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
        INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
        INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
        INFO = -8
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
        INFO = -13
      ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
        INFO = -15
      END IF
      IF( INFO.NE.0 ) THEN
        CALL XERBLA( 'DGERFSX', -INFO )
        RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) THEN
         RCOND = 1.0D+0
         DO J = 1, NRHS
            BERR( J ) = 0.0D+0
            IF ( N_ERR_BNDS .GE. 1 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            ELSE IF ( N_ERR_BNDS .GE. 2 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I) = 0.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 0.0D+0
            ELSE IF ( N_ERR_BNDS .GE. 3 ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 1.0D+0
            END IF
         END DO
         RETURN
      END IF
*
*     Default to failure.
*
      RCOND = 0.0D+0
      DO J = 1, NRHS
         BERR( J ) = 1.0D+0
         IF ( N_ERR_BNDS .GE. 1 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
         ELSE IF ( N_ERR_BNDS .GE. 2 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
         ELSE IF ( N_ERR_BNDS .GE. 3 ) THEN
            ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = 0.0D+0
            ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = 0.0D+0
         END IF
      END DO
*
*     Compute the norm of A and the reciprocal of the condition
*     number of A.
*
      IF( NOTRAN ) THEN
         NORM = 'I'
      ELSE
         NORM = '1'
      END IF
      ANORM = DLANGE( NORM, N, N, A, LDA, WORK )
      CALL DGECON( NORM, N, AF, LDAF, ANORM, RCOND, WORK, IWORK, INFO )
*
*     Perform refinement on each right-hand side
*
      IF ( REF_TYPE .NE. 0 ) THEN

         PREC_TYPE = ILAPREC( 'E' )

         IF ( NOTRAN ) THEN
            CALL DLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N,
     $           NRHS, A, LDA, AF, LDAF, IPIV, COLEQU, C, B,
     $           LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM,
     $           ERR_BNDS_COMP, WORK(N+1), WORK(1), WORK(2*N+1),
     $           WORK(1), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH,
     $           IGNORE_CWISE, INFO )
         ELSE
            CALL DLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE,  N,
     $           NRHS, A, LDA, AF, LDAF, IPIV, ROWEQU, R, B,
     $           LDB, X, LDX, BERR, N_NORMS, ERR_BNDS_NORM,
     $           ERR_BNDS_COMP, WORK(N+1), WORK(1), WORK(2*N+1),
     $           WORK(1), RCOND, ITHRESH, RTHRESH, UNSTABLE_THRESH,
     $           IGNORE_CWISE, INFO )
         END IF
      END IF

      ERR_LBND = MAX( 10.0D+0, SQRT( DBLE( N ) ) ) * DLAMCH( 'Epsilon' )
      IF ( N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 1 ) THEN
*
*     Compute scaled normwise condition number cond(A*C).
*
         IF ( COLEQU .AND. NOTRAN ) THEN
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV,
     $           -1, C, INFO, WORK, IWORK )
         ELSE IF ( ROWEQU .AND. .NOT. NOTRAN ) THEN
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV,
     $           -1, R, INFO, WORK, IWORK )
         ELSE
            RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF, IPIV,
     $           0, R, INFO, WORK, IWORK )
         END IF
         DO J = 1, NRHS
*
*     Cap the error at 1.0.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I
     $           .AND. ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .GT. 1.0D+0 )
     $           ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
*
*     Threshold the error (see LAWN).
*
            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 0.0D+0
               IF ( INFO .LE. N ) INFO = N + J
            ELSE IF ( ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) .LT. ERR_LBND )
     $     THEN
               ERR_BNDS_NORM( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_NORM( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
*
*     Save the condition number.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_NORM( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF
         END DO
      END IF

      IF ( N_ERR_BNDS .GE. 1 .AND. N_NORMS .GE. 2 ) THEN
*
*     Compute componentwise condition number cond(A*diag(Y(:,J))) for
*     each right-hand side using the current solution as an estimate of
*     the true solution.  If the componentwise error estimate is too
*     large, then the solution is a lousy estimate of truth and the
*     estimated RCOND may be too optimistic.  To avoid misleading users,
*     the inverse condition number is set to 0.0 when the estimated
*     cwise error is at least CWISE_WRONG.
*
         CWISE_WRONG = SQRT( DLAMCH( 'Epsilon' ) )
         DO J = 1, NRHS
            IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .LT. CWISE_WRONG )
     $           THEN
               RCOND_TMP = DLA_GERCOND( TRANS, N, A, LDA, AF, LDAF,
     $              IPIV, 1, X(1,J), INFO, WORK, IWORK )
            ELSE
               RCOND_TMP = 0.0D+0
            END IF
*
*     Cap the error at 1.0.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_ERR_I
     $           .AND. ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) .GT. 1.0D+0 )
     $           ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
*
*     Threshold the error (see LAWN).
*
            IF ( RCOND_TMP .LT. ILLRCOND_THRESH ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = 1.0D+0
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 0.0D+0
               IF ( PARAMS( LA_LINRX_CWISE_I ) .EQ. 1.0D+0
     $              .AND. INFO.LT.N + J ) INFO = N + J
            ELSE IF ( ERR_BNDS_COMP( J, LA_LINRX_ERR_I )
     $              .LT. ERR_LBND ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_ERR_I ) = ERR_LBND
               ERR_BNDS_COMP( J, LA_LINRX_TRUST_I ) = 1.0D+0
            END IF
*
*     Save the condition number.
*
            IF ( N_ERR_BNDS .GE. LA_LINRX_RCOND_I ) THEN
               ERR_BNDS_COMP( J, LA_LINRX_RCOND_I ) = RCOND_TMP
            END IF
         END DO
      END IF
*
      RETURN
*
*     End of DGERFSX
*
      END
