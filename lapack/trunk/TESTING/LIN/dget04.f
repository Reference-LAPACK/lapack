      SUBROUTINE DGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            LDX, LDXACT, N, NRHS
      DOUBLE PRECISION   RCOND, RESID
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( LDX, * ), XACT( LDXACT, * )
*     ..
*
*  Purpose
*  =======
*
*  DGET04 computes the difference between a computed solution and the
*  true solution to a system of linear equations.
*
*  RESID =  ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),
*  where RCOND is the reciprocal of the condition number and EPS is the
*  machine epsilon.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of the matrices X and XACT.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of columns of the matrices X and XACT.  NRHS >= 0.
*
*  X       (input) DOUBLE PRECISION array, dimension (LDX,NRHS)
*          The computed solution vectors.  Each vector is stored as a
*          column of the matrix X.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X.  LDX >= max(1,N).
*
*  XACT    (input) DOUBLE PRECISION array, dimension( LDX, NRHS )
*          The exact solution vectors.  Each vector is stored as a
*          column of the matrix XACT.
*
*  LDXACT  (input) INTEGER
*          The leading dimension of the array XACT.  LDXACT >= max(1,N).
*
*  RCOND   (input) DOUBLE PRECISION
*          The reciprocal of the condition number of the coefficient
*          matrix in the system of equations.
*
*  RESID   (output) DOUBLE PRECISION
*          The maximum over the NRHS solution vectors of
*          ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS )
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IX, J
      DOUBLE PRECISION   DIFFNM, EPS, XNORM
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           IDAMAX, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Quick exit if N = 0 or NRHS = 0.
*
      IF( N.LE.0 .OR. NRHS.LE.0 ) THEN
         RESID = ZERO
         RETURN
      END IF
*
*     Exit with RESID = 1/EPS if RCOND is invalid.
*
      EPS = DLAMCH( 'Epsilon' )
      IF( RCOND.LT.ZERO ) THEN
         RESID = 1.0D0 / EPS
         RETURN
      END IF
*
*     Compute the maximum of
*        norm(X - XACT) / ( norm(XACT) * EPS )
*     over all the vectors X and XACT .
*
      RESID = ZERO
      DO 20 J = 1, NRHS
         IX = IDAMAX( N, XACT( 1, J ), 1 )
         XNORM = ABS( XACT( IX, J ) )
         DIFFNM = ZERO
         DO 10 I = 1, N
            DIFFNM = MAX( DIFFNM, ABS( X( I, J )-XACT( I, J ) ) )
   10    CONTINUE
         IF( XNORM.LE.ZERO ) THEN
            IF( DIFFNM.GT.ZERO )
     $         RESID = 1.0D0 / EPS
         ELSE
            RESID = MAX( RESID, ( DIFFNM / XNORM )*RCOND )
         END IF
   20 CONTINUE
      IF( RESID*EPS.LT.1.0D0 )
     $   RESID = RESID / EPS
*
      RETURN
*
*     End of DGET04
*
      END
