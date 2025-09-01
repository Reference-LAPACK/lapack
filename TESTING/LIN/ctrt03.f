*> \brief \b CTRT03
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CTRT03( UPLO, TRANS, DIAG, N, NRHS, A, LDA, SCALE,
*                          CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, TRANS, UPLO
*       INTEGER            LDA, LDB, LDX, N, NRHS
*       REAL               RESID, SCALE, TSCAL
*       ..
*       .. Array Arguments ..
*       REAL               CNORM( * )
*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ),
*      $                   X( LDX, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTRT03 computes the residual for the solution to a scaled triangular
*> system of equations A*x = s*b,  A**T *x = s*b,  or  A**H *x = s*b.
*> Here A is a triangular matrix, A**T denotes the transpose of A, A**H
*> denotes the conjugate transpose of A, s is a scalar, and x and b are
*> N by NRHS matrices.  The test ratio is the maximum over the number of
*> right hand sides of
*>    norm(s*b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
*> where op(A) denotes A, A**T, or A**H, and EPS is the machine epsilon.
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
*>          = 'N':  A *x = s*b     (No transpose)
*>          = 'T':  A**T *x = s*b  (Transpose)
*>          = 'C':  A**H *x = s*b  (Conjugate transpose)
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
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          The triangular matrix A.  If UPLO = 'U', the leading n by n
*>          upper triangular part of the array A contains the upper
*>          triangular matrix, and the strictly lower triangular part of
*>          A is not referenced.  If UPLO = 'L', the leading n by n lower
*>          triangular part of the array A contains the lower triangular
*>          matrix, and the strictly upper triangular part of A is not
*>          referenced.  If DIAG = 'U', the diagonal elements of A are
*>          also not referenced and are assumed to be 1.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] SCALE
*> \verbatim
*>          SCALE is REAL
*>          The scaling factor s used in solving the triangular system.
*> \endverbatim
*>
*> \param[in] CNORM
*> \verbatim
*>          CNORM is REAL array, dimension (N)
*>          The 1-norms of the columns of A, not counting the diagonal.
*> \endverbatim
*>
*> \param[in] TSCAL
*> \verbatim
*>          TSCAL is REAL
*>          The scaling factor used in computing the 1-norms in CNORM.
*>          CNORM actually contains the column norms of TSCAL*A.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX array, dimension (LDX,NRHS)
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
*>          B is COMPLEX array, dimension (LDB,NRHS)
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
*>          WORK is COMPLEX array, dimension (N)
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
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
*> \ingroup complex_lin
*
*  =====================================================================
      SUBROUTINE CTRT03( UPLO, TRANS, DIAG, N, NRHS, A, LDA, SCALE,
     $                   CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            LDA, LDB, LDX, N, NRHS
      REAL               RESID, SCALE, TSCAL
*     ..
*     .. Array Arguments ..
      REAL               CNORM( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ),
     $                   X( LDX, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX, J
      REAL               EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ICAMAX
      REAL               SLAMCH
      EXTERNAL           LSAME, ICAMAX, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CAXPY, CCOPY, CSSCAL, CTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, MAX, REAL
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
      EPS = SLAMCH( 'Epsilon' )
      SMLNUM = SLAMCH( 'Safe minimum' )
*
*     Compute the norm of the triangular matrix A using the column
*     norms already computed by CLATRS.
*
      TNORM = ZERO
      IF( LSAME( DIAG, 'N' ) ) THEN
         DO 10 J = 1, N
            TNORM = MAX( TNORM, TSCAL*ABS( A( J, J ) )+CNORM( J ) )
   10    CONTINUE
      ELSE
         DO 20 J = 1, N
            TNORM = MAX( TNORM, TSCAL+CNORM( J ) )
   20    CONTINUE
      END IF
*
*     Compute the maximum over the number of right hand sides of
*        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
*
      RESID = ZERO
      DO 30 J = 1, NRHS
         CALL CCOPY( N, X( 1, J ), 1, WORK, 1 )
         IX = ICAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / REAL( N )
         CALL CSSCAL( N, XSCAL, WORK, 1 )
         CALL CTRMV( UPLO, TRANS, DIAG, N, A, LDA, WORK, 1 )
         CALL CAXPY( N, CMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 )
         IX = ICAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = ICAMAX( N, X( 1, J ), 1 )
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
   30 CONTINUE
*
      RETURN
*
*     End of CTRT03
*
      END
