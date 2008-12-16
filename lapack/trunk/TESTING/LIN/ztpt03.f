      SUBROUTINE ZTPT03( UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM,
     $                   TSCAL, X, LDX, B, LDB, WORK, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            LDB, LDX, N, NRHS
      DOUBLE PRECISION   RESID, SCALE, TSCAL
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   CNORM( * )
      COMPLEX*16         AP( * ), B( LDB, * ), WORK( * ), X( LDX, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTPT03 computes the residual for the solution to a scaled triangular
*  system of equations A*x = s*b,  A**T *x = s*b,  or  A**H *x = s*b,
*  when the triangular matrix A is stored in packed format.  Here A**T
*  denotes the transpose of A, A**H denotes the conjugate transpose of
*  A, s is a scalar, and x and b are N by NRHS matrices.  The test ratio
*  is the maximum over the number of right hand sides of
*     norm(s*b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
*  where op(A) denotes A, A**T, or A**H, and EPS is the machine epsilon.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  TRANS   (input) CHARACTER*1
*          Specifies the operation applied to A.
*          = 'N':  A *x = s*b     (No transpose)
*          = 'T':  A**T *x = s*b  (Transpose)
*          = 'C':  A**H *x = s*b  (Conjugate transpose)
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrices X and B.  NRHS >= 0.
*
*  AP      (input) COMPLEX*16 array, dimension (N*(N+1)/2)
*          The upper or lower triangular matrix A, packed columnwise in
*          a linear array.  The j-th column of A is stored in the array
*          AP as follows:
*          if UPLO = 'U', AP((j-1)*j/2 + i) = A(i,j) for 1<=i<=j;
*          if UPLO = 'L',
*             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = A(i,j) for j<=i<=n.
*
*  SCALE   (input) DOUBLE PRECISION
*          The scaling factor s used in solving the triangular system.
*
*  CNORM   (input) DOUBLE PRECISION array, dimension (N)
*          The 1-norms of the columns of A, not counting the diagonal.
*
*  TSCAL   (input) DOUBLE PRECISION
*          The scaling factor used in computing the 1-norms in CNORM.
*          CNORM actually contains the column norms of TSCAL*A.
*
*  X       (input) COMPLEX*16 array, dimension (LDX,NRHS)
*          The computed solution vectors for the system of linear
*          equations.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  B       (input) COMPLEX*16 array, dimension (LDB,NRHS)
*          The right hand side vectors for the system of linear
*          equations.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  RESID   (output) DOUBLE PRECISION
*          The maximum over the number of right hand sides of
*          norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX, J, JJ
      DOUBLE PRECISION   EPS, ERR, SMLNUM, TNORM, XNORM, XSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IZAMAX, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZCOPY, ZDSCAL, ZTPMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, MAX
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
*
*     Compute the norm of the triangular matrix A using the column
*     norms already computed by ZLATPS.
*
      TNORM = 0.D0
      IF( LSAME( DIAG, 'N' ) ) THEN
         IF( LSAME( UPLO, 'U' ) ) THEN
            JJ = 1
            DO 10 J = 1, N
               TNORM = MAX( TNORM, TSCAL*ABS( AP( JJ ) )+CNORM( J ) )
               JJ = JJ + J
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
*        norm(op(A)*x - s*b) / ( norm(A) * norm(x) * EPS ).
*
      RESID = ZERO
      DO 40 J = 1, NRHS
         CALL ZCOPY( N, X( 1, J ), 1, WORK, 1 )
         IX = IZAMAX( N, WORK, 1 )
         XNORM = MAX( ONE, ABS( X( IX, J ) ) )
         XSCAL = ( ONE / XNORM ) / DBLE( N )
         CALL ZDSCAL( N, XSCAL, WORK, 1 )
         CALL ZTPMV( UPLO, TRANS, DIAG, N, AP, WORK, 1 )
         CALL ZAXPY( N, DCMPLX( -SCALE*XSCAL ), B( 1, J ), 1, WORK, 1 )
         IX = IZAMAX( N, WORK, 1 )
         ERR = TSCAL*ABS( WORK( IX ) )
         IX = IZAMAX( N, X( 1, J ), 1 )
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
*     End of ZTPT03
*
      END
