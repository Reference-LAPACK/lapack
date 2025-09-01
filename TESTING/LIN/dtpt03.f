*> \brief \b DTPT03
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTPT03( UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM,
*                          TSCAL, X, LDX, B, LDB, WORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, TRANS, UPLO
*       INTEGER            LDB, LDX, N, NRHS
*       DOUBLE PRECISION   RESID, SCALE, TSCAL
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   AP( * ), B( LDB, * ), CNORM( * ), WORK( * ),
*      $                   X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTPT03 computes the residual for the solution to a scaled triangular
*> system of equations A*x = s*b  or  A'*x = s*b  when the triangular
*> matrix A is stored in packed format.  Here A' is the transpose of A,
*> s is a scalar, and x and b are N by NRHS matrices.  The test ratio is
*> the maximum over the number of right hand sides of
*>    norm(s*b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
*> where op(A) denotes A or A' and EPS is the machine epsilon.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the matrix A is upper or lower triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the operation applied to A.
*>          = 'N':  A *x = s*b  (No transpose)
*>          = 'T':  A'*x = s*b  (Transpose)
*>          = 'C':  A'*x = s*b  (Conjugate transpose = Transpose)
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the matrix A is unit triangular.
*>          = 'N':  Non-unit triangular
*>          = 'U':  Unit triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrices X and B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] AP
*> \verbatim
*>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
*>          The upper or lower triangular matrix A, packed columnwise in
*>          a linear array.  The j-th column of A is stored in the array
*>          AP as follows:
*>          if UPLO = 'U', AP((j-1)*j/2 + i) = A(i,j) for 1<=i<=j;
*>          if UPLO = 'L',
*>             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = A(i,j) for j<=i<=n.
*> \endverbatim
*>
*> \param[in] SCALE
*> \verbatim
*>          SCALE is DOUBLE PRECISION
*>          The scaling factor s used in solving the triangular system.
*> \endverbatim
*>
*> \param[in] CNORM
*> \verbatim
*>          CNORM is DOUBLE PRECISION array, dimension (N)
*>          The 1-norms of the columns of A, not counting the diagonal.
*> \endverbatim
*>
*> \param[in] TSCAL
*> \verbatim
*>          TSCAL is DOUBLE PRECISION
*>          The scaling factor used in computing the 1-norms in CNORM.
*>          CNORM actually contains the column norms of TSCAL*A.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
*>          The computed solution vectors for the system of linear
*>          equations.
*> \endverbatim
*>
*> \param[in] LDX
*> \verbatim
*>          LDX is INTEGER
*>          The leading dimension of the array X.  LDX >= max(1,N).
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          The right hand side vectors for the system of linear
*>          equations.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is DOUBLE PRECISION
*>          The maximum over the number of right hand sides of
*>          norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DTPT03( UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM,
     $                   TSCAL, X, LDX, B, LDB, WORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID, SCALE, TSCAL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   AP( * ), B( LDB, * ), CNORM( * ), WORK( * ),
     $                   X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX, J, JJ
      DOUBLE PRECISION   BIGNUM, EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DCOPY, DSCAL, DTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
      EPS = DLAMCH( 'Epsilon' )
      SMLNUM = DLAMCH( 'Safe minimum' )
      BIGNUM = ONE / SMLNUM
*
*     Compute the norm of the triangular matrix A using the column
*     norms already computed by DLATPS.
*
      TNORM = ZERO
      IF( LSAME( DIAG, 'N' ) ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            JJ = 1
            DO 10 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + J + 1
   10       CONTINUE
         ELSE
            JJ = 1
            DO 20 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + N - J + 1
   20       CONTINUE
         END IF
      ELSE
         DO 30 J = 1, N
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
   30    CONTINUE
      END IF
*
*     Compute the maximum over the number of right hand sides of
*        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
*
      RESID = ZERO
      DO 40 J = 1, NRHS
         CALL DCOPY( N, X( 1, J ), 1, WORK, 1 )
         IX = IDAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / DBLE( N )
         CALL DSCAL( N, XSCAL, WORK, 1 )
         CALL DTPMV( UPLO, TRANS, DIAG, N, AP, WORK, 1 )
         CALL DAXPY( N, -SCALE*XSCAL, B( 1, J ), 1, WORK, 1 )
         IX = IDAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = IDAMAX( N, X( 1, J ), 1 )
         XNORM = ABS( X( IX, J ) )
         IF( ERR*SMLNUM.LE.XNORM ) THEN
            IF( XNORM.GT.ZERO )
     $         ERR = ERR / XNORM
         ELSE
            IF( ERR.GT.ZERO )
     $         ERR = ONE / EPS
         END IF
         IF( ERR*SMLNUM.LE.TNORM ) THEN
            IF( TNORM.GT.ZERO )
     $         ERR = ERR / TNORM
         ELSE
            IF( ERR.GT.ZERO )
     $         ERR = ONE / EPS
         END IF
         RESID = MAX( RESID, ERR )
   40 CONTINUE
*
      RETURN
*
*     End of DTPT03
*
      END
