*> \brief <b> CPOSVXX computes the solution to system of linear equations A * X = B for PO matrices</b>
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CPOSVXX + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cposvxx.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cposvxx.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cposvxx.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, EQUED,
*                           S, B, LDB, X, LDX, RCOND, RPVGRW, BERR,
*                           N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP,
*                           NPARAMS, PARAMS, WORK, RWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          EQUED, FACT, UPLO
*       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,
*      $                   N_ERR_BNDS
*       REAL               RCOND, RPVGRW
*       ..
*       .. Array Arguments ..
*       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
*      $                   WORK( * ), X( LDX, * )
*       REAL               S( * ), PARAMS( * ), BERR( * ), RWORK( * ),
*      $                   ERR_BNDS_NORM( NRHS, * ),
*      $                   ERR_BNDS_COMP( NRHS, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    CPOSVXX uses the Cholesky factorization A = U**T*U or A = L*L**T
*>    to compute the solution to a complex system of linear equations
*>    A * X = B, where A is an N-by-N Hermitian positive definite matrix
*>    and X and B are N-by-NRHS matrices.
*>
*>    If requested, both normwise and maximum componentwise error bounds
*>    are returned. CPOSVXX will return a solution with a tiny
*>    guaranteed error (O(eps) where eps is the working machine
*>    precision) unless the matrix is very ill-conditioned, in which
*>    case a warning is returned. Relevant condition numbers also are
*>    calculated and returned.
*>
*>    CPOSVXX accepts user-provided factorizations and equilibration
*>    factors; see the definitions of the FACT and EQUED options.
*>    Solving with refinement and using a factorization from a previous
*>    CPOSVXX call will also produce a solution with either O(eps)
*>    errors or warnings, but we cannot make that claim for general
*>    user-provided factorizations and equilibration factors if they
*>    differ from what CPOSVXX would itself produce.
*> \endverbatim
*
*> \par Description:
*  =================
*>
*> \verbatim
*>
*>    The following steps are performed:
*>
*>    1. If FACT = 'E', real scaling factors are computed to equilibrate
*>    the system:
*>
*>      diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B
*>
*>    Whether or not the system will be equilibrated depends on the
*>    scaling of the matrix A, but if equilibration is used, A is
*>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.
*>
*>    2. If FACT = 'N' or 'E', the Cholesky decomposition is used to
*>    factor the matrix A (after equilibration if FACT = 'E') as
*>       A = U**T* U,  if UPLO = 'U', or
*>       A = L * L**T,  if UPLO = 'L',
*>    where U is an upper triangular matrix and L is a lower triangular
*>    matrix.
*>
*>    3. If the leading principal minor of order i is not positive,
*>    then the routine returns with INFO = i. Otherwise, the factored
*>    form of A is used to estimate the condition number of the matrix
*>    A (see argument RCOND).  If the reciprocal of the condition number
*>    is less than machine precision, the routine still goes on to solve
*>    for X and compute error bounds as described below.
*>
*>    4. The system of equations is solved for X using the factored form
*>    of A.
*>
*>    5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),
*>    the routine will use iterative refinement to try to get a small
*>    error and error bounds.  Refinement calculates the residual to at
*>    least twice the working precision.
*>
*>    6. If equilibration was used, the matrix X is premultiplied by
*>    diag(S) so that it solves the original system before
*>    equilibration.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \verbatim
*>     Some optional parameters are bundled in the PARAMS array.  These
*>     settings determine how refinement is performed, but often the
*>     defaults are acceptable.  If the defaults are acceptable, users
*>     can pass NPARAMS = 0 which prevents the source code from accessing
*>     the PARAMS argument.
*> \endverbatim
*>
*> \param[in] FACT
*> \verbatim
*>          FACT is CHARACTER*1
*>     Specifies whether or not the factored form of the matrix A is
*>     supplied on entry, and if not, whether the matrix A should be
*>     equilibrated before it is factored.
*>       = 'F':  On entry, AF contains the factored form of A.
*>               If EQUED is not 'N', the matrix A has been
*>               equilibrated with scaling factors given by S.
*>               A and AF are not modified.
*>       = 'N':  The matrix A will be copied to AF and factored.
*>       = 'E':  The matrix A will be equilibrated if necessary, then
*>               copied to AF and factored.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>       = 'U':  Upper triangle of A is stored;
*>       = 'L':  Lower triangle of A is stored.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>     The number of linear equations, i.e., the order of the
*>     matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>     The number of right hand sides, i.e., the number of columns
*>     of the matrices B and X.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>     On entry, the Hermitian matrix A, except if FACT = 'F' and EQUED =
*>     'Y', then A must contain the equilibrated matrix
*>     diag(S)*A*diag(S).  If UPLO = 'U', the leading N-by-N upper
*>     triangular part of A contains the upper triangular part of the
*>     matrix A, and the strictly lower triangular part of A is not
*>     referenced.  If UPLO = 'L', the leading N-by-N lower triangular
*>     part of A contains the lower triangular part of the matrix A, and
*>     the strictly upper triangular part of A is not referenced.  A is
*>     not modified if FACT = 'F' or 'N', or if FACT = 'E' and EQUED =
*>     'N' on exit.
*>
*>     On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
*>     diag(S)*A*diag(S).
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>     The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] AF
*> \verbatim
*>          AF is COMPLEX array, dimension (LDAF,N)
*>     If FACT = 'F', then AF is an input argument and on entry
*>     contains the triangular factor U or L from the Cholesky
*>     factorization A = U**T*U or A = L*L**T, in the same storage
*>     format as A.  If EQUED .ne. 'N', then AF is the factored
*>     form of the equilibrated matrix diag(S)*A*diag(S).
*>
*>     If FACT = 'N', then AF is an output argument and on exit
*>     returns the triangular factor U or L from the Cholesky
*>     factorization A = U**T*U or A = L*L**T of the original
*>     matrix A.
*>
*>     If FACT = 'E', then AF is an output argument and on exit
*>     returns the triangular factor U or L from the Cholesky
*>     factorization A = U**T*U or A = L*L**T of the equilibrated
*>     matrix A (see the description of A for the form of the
*>     equilibrated matrix).
*> \endverbatim
*>
*> \param[in] LDAF
*> \verbatim
*>          LDAF is INTEGER
*>     The leading dimension of the array AF.  LDAF >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] EQUED
*> \verbatim
*>          EQUED is CHARACTER*1
*>     Specifies the form of equilibration that was done.
*>       = 'N':  No equilibration (always true if FACT = 'N').
*>       = 'Y':  Both row and column equilibration, i.e., A has been
*>               replaced by diag(S) * A * diag(S).
*>     EQUED is an input argument if FACT = 'F'; otherwise, it is an
*>     output argument.
*> \endverbatim
*>
*> \param[in,out] S
*> \verbatim
*>          S is REAL array, dimension (N)
*>     The row scale factors for A.  If EQUED = 'Y', A is multiplied on
*>     the left and right by diag(S).  S is an input argument if FACT =
*>     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED
*>     = 'Y', each element of S must be positive.  If S is output, each
*>     element of S is a power of the radix. If S is input, each element
*>     of S should be a power of the radix to ensure a reliable solution
*>     and error estimates. Scaling by powers of the radix does not cause
*>     rounding errors unless the result underflows or overflows.
*>     Rounding errors during scaling lead to refining with a matrix that
*>     is not equivalent to the input matrix, producing error estimates
*>     that may not be reliable.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,NRHS)
*>     On entry, the N-by-NRHS right hand side matrix B.
*>     On exit,
*>     if EQUED = 'N', B is not modified;
*>     if EQUED = 'Y', B is overwritten by diag(S)*B;
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>     The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX array, dimension (LDX,NRHS)
*>     If INFO = 0, the N-by-NRHS solution matrix X to the original
*>     system of equations.  Note that A and B are modified on exit if
*>     EQUED .ne. 'N', and the solution to the equilibrated system is
*>     inv(diag(S))*X.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>     The leading dimension of the array X.  LDX >= max(1,N).
*> \endverbatim
*>
*> \param[out] RCOND
*> \verbatim
*>          RCOND is REAL
*>     Reciprocal scaled condition number.  This is an estimate of the
*>     reciprocal Skeel condition number of the matrix A after
*>     equilibration (if done).  If this is less than the machine
*>     precision (in particular, if it is zero), the matrix is singular
*>     to working precision.  Note that the error may still be small even
*>     if this number is very small and the matrix appears ill-
*>     conditioned.
*> \endverbatim
*>
*> \param[out] RPVGRW
*> \verbatim
*>          RPVGRW is REAL
*>     Reciprocal pivot growth.  On exit, this contains the reciprocal
*>     pivot growth factor norm(A)/norm(U). The "max absolute element"
*>     norm is used.  If this is much less than 1, then the stability of
*>     the LU factorization of the (equilibrated) matrix A could be poor.
*>     This also means that the solution X, estimated condition numbers,
*>     and error bounds could be unreliable. If factorization fails with
*>     0<INFO<=N, then this contains the reciprocal pivot growth factor
*>     for the leading INFO columns of A.
*> \endverbatim
*>
*> \param[out] BERR
*> \verbatim
*>          BERR is REAL array, dimension (NRHS)
*>     Componentwise relative backward error.  This is the
*>     componentwise relative backward error of each solution vector X(j)
*>     (i.e., the smallest relative change in any element of A or B that
*>     makes X(j) an exact solution).
*> \endverbatim
*>
*> \param[in] N_ERR_BNDS
*> \verbatim
*>          N_ERR_BNDS is INTEGER
*>     Number of error bounds to return for each right hand side
*>     and each type (normwise or componentwise).  See ERR_BNDS_NORM and
*>     ERR_BNDS_COMP below.
*> \endverbatim
*>
*> \param[out] ERR_BNDS_NORM
*> \verbatim
*>          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)
*>     For each right-hand side, this array contains information about
*>     various error bounds and condition numbers corresponding to the
*>     normwise relative error, which is defined as follows:
*>
*>     Normwise relative error in the ith solution vector:
*>             max_j (abs(XTRUE(j,i) - X(j,i)))
*>            ------------------------------
*>                  max_j abs(X(j,i))
*>
*>     The array is indexed by the type of error information as described
*>     below. There currently are up to three pieces of information
*>     returned.
*>
*>     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
*>     right-hand side.
*>
*>     The second index in ERR_BNDS_NORM(:,err) contains the following
*>     three fields:
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the
*>              reciprocal condition number is less than the threshold
*>              sqrt(n) * slamch('Epsilon').
*>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,
*>              almost certainly within a factor of 10 of the true error
*>              so long as the next entry is greater than the threshold
*>              sqrt(n) * slamch('Epsilon'). This error bound should only
*>              be trusted if the previous boolean is true.
*>
*>     err = 3  Reciprocal condition number: Estimated normwise
*>              reciprocal condition number.  Compared with the threshold
*>              sqrt(n) * slamch('Epsilon') to determine if the error
*>              estimate is "guaranteed". These reciprocal condition
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
*>              appropriately scaled matrix Z.
*>              Let Z = S*A, where S scales each row by a power of the
*>              radix so all absolute row sums of Z are approximately 1.
*>
*>     See Lapack Working Note 165 for further details and extra
*>     cautions.
*> \endverbatim
*>
*> \param[out] ERR_BNDS_COMP
*> \verbatim
*>          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS)
*>     For each right-hand side, this array contains information about
*>     various error bounds and condition numbers corresponding to the
*>     componentwise relative error, which is defined as follows:
*>
*>     Componentwise relative error in the ith solution vector:
*>                    abs(XTRUE(j,i) - X(j,i))
*>             max_j ----------------------
*>                         abs(X(j,i))
*>
*>     The array is indexed by the right-hand side i (on which the
*>     componentwise relative error depends), and the type of error
*>     information as described below. There currently are up to three
*>     pieces of information returned for each right-hand side. If
*>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then
*>     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
*>     the first (:,N_ERR_BNDS) entries are returned.
*>
*>     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
*>     right-hand side.
*>
*>     The second index in ERR_BNDS_COMP(:,err) contains the following
*>     three fields:
*>     err = 1 "Trust/don't trust" boolean. Trust the answer if the
*>              reciprocal condition number is less than the threshold
*>              sqrt(n) * slamch('Epsilon').
*>
*>     err = 2 "Guaranteed" error bound: The estimated forward error,
*>              almost certainly within a factor of 10 of the true error
*>              so long as the next entry is greater than the threshold
*>              sqrt(n) * slamch('Epsilon'). This error bound should only
*>              be trusted if the previous boolean is true.
*>
*>     err = 3  Reciprocal condition number: Estimated componentwise
*>              reciprocal condition number.  Compared with the threshold
*>              sqrt(n) * slamch('Epsilon') to determine if the error
*>              estimate is "guaranteed". These reciprocal condition
*>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
*>              appropriately scaled matrix Z.
*>              Let Z = S*(A*diag(x)), where x is the solution for the
*>              current right-hand side and S scales each row of
*>              A*diag(x) by a power of the radix so all absolute row
*>              sums of Z are approximately 1.
*>
*>     See Lapack Working Note 165 for further details and extra
*>     cautions.
*> \endverbatim
*>
*> \param[in] NPARAMS
*> \verbatim
*>          NPARAMS is INTEGER
*>     Specifies the number of parameters set in PARAMS.  If <= 0, the
*>     PARAMS array is never referenced and default values are used.
*> \endverbatim
*>
*> \param[in,out] PARAMS
*> \verbatim
*>          PARAMS is REAL array, dimension NPARAMS
*>     Specifies algorithm parameters.  If an entry is < 0.0, then
*>     that entry will be filled with default value used for that
*>     parameter.  Only positions up to NPARAMS are accessed; defaults
*>     are used for higher-numbered parameters.
*>
*>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
*>            refinement or not.
*>         Default: 1.0
*>            = 0.0:  No refinement is performed, and no error bounds are
*>                    computed.
*>            = 1.0:  Use the double-precision refinement algorithm,
*>                    possibly with doubled-single computations if the
*>                    compilation environment does not support DOUBLE
*>                    PRECISION.
*>              (other values are reserved for future use)
*>
*>       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
*>            computations allowed for refinement.
*>         Default: 10
*>         Aggressive: Set to 100 to permit convergence using approximate
*>                     factorizations or factorizations other than LU. If
*>                     the factorization uses a technique other than
*>                     Gaussian elimination, the guarantees in
*>                     err_bnds_norm and err_bnds_comp may no longer be
*>                     trustworthy.
*>
*>       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
*>            will attempt to find a solution with small componentwise
*>            relative error in the double-precision algorithm.  Positive
*>            is true, 0.0 is false.
*>         Default: 1.0 (attempt componentwise convergence)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (2*N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>       = 0:  Successful exit. The solution to every right-hand side is
*>         guaranteed.
*>       < 0:  If INFO = -i, the i-th argument had an illegal value
*>       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
*>         has been completed, but the factor U is exactly singular, so
*>         the solution and error bounds could not be computed. RCOND = 0
*>         is returned.
*>       = N+J: The solution corresponding to the Jth right-hand side is
*>         not guaranteed. The solutions corresponding to other right-
*>         hand sides K with K > J may not be guaranteed as well, but
*>         only the first such right-hand side is reported. If a small
*>         componentwise error is not requested (PARAMS(3) = 0.0) then
*>         the Jth right-hand side is the first with a normwise error
*>         bound that is not guaranteed (the smallest J such
*>         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
*>         the Jth right-hand side is the first with either a normwise or
*>         componentwise error bound that is not guaranteed (the smallest
*>         J such that either ERR_BNDS_NORM(J,1) = 0.0 or
*>         ERR_BNDS_COMP(J,1) = 0.0). See the definition of
*>         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
*>         about all of the right-hand sides check ERR_BNDS_NORM or
*>         ERR_BNDS_COMP.
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
*> \ingroup posvxx
*
*  =====================================================================
      SUBROUTINE CPOSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF,
     $                    EQUED,
     $                    S, B, LDB, X, LDX, RCOND, RPVGRW, BERR,
     $                    N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP,
     $                    NPARAMS, PARAMS, WORK, RWORK, INFO )
*
*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED, FACT, UPLO
      INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,
     $                   N_ERR_BNDS
      REAL               RCOND, RPVGRW
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
     $                   WORK( * ), X( LDX, * )
      REAL               S( * ), PARAMS( * ), BERR( * ), RWORK( * ),
     $                   ERR_BNDS_NORM( NRHS, * ),
     $                   ERR_BNDS_COMP( NRHS, * )
*     ..
*
*  ==================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      INTEGER            FINAL_NRM_ERR_I, FINAL_CMP_ERR_I, BERR_I
      INTEGER            RCOND_I, NRM_RCOND_I, NRM_ERR_I, CMP_RCOND_I
      INTEGER            CMP_ERR_I, PIV_GROWTH_I
      PARAMETER          ( FINAL_NRM_ERR_I = 1, FINAL_CMP_ERR_I = 2,
     $                   BERR_I = 3 )
      PARAMETER          ( RCOND_I = 4, NRM_RCOND_I = 5, NRM_ERR_I = 6 )
      PARAMETER          ( CMP_RCOND_I = 7, CMP_ERR_I = 8,
     $                   PIV_GROWTH_I = 9 )
*     ..
*     .. Local Scalars ..
      LOGICAL            EQUIL, NOFACT, RCEQU
      INTEGER            INFEQU, J
      REAL               AMAX, BIGNUM, SMIN, SMAX, SCOND, SMLNUM
*     ..
*     .. External Functions ..
      EXTERNAL           LSAME, SLAMCH, CLA_PORPVGRW
      LOGICAL            LSAME
      REAL               SLAMCH, CLA_PORPVGRW
*     ..
*     .. External Subroutines ..
      EXTERNAL           CPOEQUB, CPOTRF, CPOTRS, CLACPY,
     $                   CLAQHE, XERBLA, CLASCL2, CPORFSX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      NOFACT = LSAME( FACT, 'N' )
      EQUIL = LSAME( FACT, 'E' )
      SMLNUM = SLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
      IF( NOFACT .OR. EQUIL ) THEN
         EQUED = 'N'
         RCEQU = .FALSE.
      ELSE
         RCEQU = LSAME( EQUED, 'Y' )
      ENDIF
*
*     Default is failure.  If an input parameter is wrong or
*     factorization fails, make everything look horrible.  Only the
*     pivot growth is set here, the rest is initialized in CPORFSX.
*
      RPVGRW = ZERO
*
*     Test the input parameters.  PARAMS is not tested until CPORFSX.
*
      IF( .NOT.NOFACT .AND. .NOT.EQUIL .AND. .NOT.
     $     LSAME( FACT, 'F' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( UPLO, 'U' ) .AND.
     $         .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDAF.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LSAME( FACT, 'F' ) .AND. .NOT.
     $        ( RCEQU .OR. LSAME( EQUED, 'N' ) ) ) THEN
         INFO = -9
      ELSE
         IF ( RCEQU ) THEN
            SMIN = BIGNUM
            SMAX = ZERO
            DO 10 J = 1, N
               SMIN = MIN( SMIN, S( J ) )
               SMAX = MAX( SMAX, S( J ) )
 10         CONTINUE
            IF( SMIN.LE.ZERO ) THEN
               INFO = -10
            ELSE IF( N.GT.0 ) THEN
               SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )
            ELSE
               SCOND = ONE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            IF( LDB.LT.MAX( 1, N ) ) THEN
               INFO = -12
            ELSE IF( LDX.LT.MAX( 1, N ) ) THEN
               INFO = -14
            END IF
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPOSVXX', -INFO )
         RETURN
      END IF
*
      IF( EQUIL ) THEN
*
*     Compute row and column scalings to equilibrate the matrix A.
*
         CALL CPOEQUB( N, A, LDA, S, SCOND, AMAX, INFEQU )
         IF( INFEQU.EQ.0 ) THEN
*
*     Equilibrate the matrix.
*
            CALL CLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )
            RCEQU = LSAME( EQUED, 'Y' )
         END IF
      END IF
*
*     Scale the right-hand side.
*
      IF( RCEQU ) CALL CLASCL2( N, NRHS, S, B, LDB )
*
      IF( NOFACT .OR. EQUIL ) THEN
*
*        Compute the Cholesky factorization of A.
*
         CALL CLACPY( UPLO, N, N, A, LDA, AF, LDAF )
         CALL CPOTRF( UPLO, N, AF, LDAF, INFO )
*
*        Return if INFO is non-zero.
*
         IF( INFO.GT.0 ) THEN
*
*           Pivot in column INFO is exactly 0
*           Compute the reciprocal pivot growth factor of the
*           leading rank-deficient INFO columns of A.
*
            RPVGRW = CLA_PORPVGRW( UPLO, N, A, LDA, AF, LDAF, RWORK )
            RETURN
         END IF
      END IF
*
*     Compute the reciprocal pivot growth factor RPVGRW.
*
      RPVGRW = CLA_PORPVGRW( UPLO, N, A, LDA, AF, LDAF, RWORK )
*
*     Compute the solution matrix X.
*
      CALL CLACPY( 'Full', N, NRHS, B, LDB, X, LDX )
      CALL CPOTRS( UPLO, N, NRHS, AF, LDAF, X, LDX, INFO )
*
*     Use iterative refinement to improve the computed solution and
*     compute error bounds and backward error estimates for it.
*
      CALL CPORFSX( UPLO, EQUED, N, NRHS, A, LDA, AF, LDAF,
     $     S, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS, ERR_BNDS_NORM,
     $     ERR_BNDS_COMP,  NPARAMS, PARAMS, WORK, RWORK, INFO )

*
*     Scale solutions.
*
      IF ( RCEQU ) THEN
         CALL CLASCL2( N, NRHS, S, X, LDX )
      END IF
*
      RETURN
*
*     End of CPOSVXX
*
      END
