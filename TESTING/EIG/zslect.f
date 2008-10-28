      LOGICAL          FUNCTION ZSLECT( Z )
*
*  -- LAPACK test routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*     .. Scalar Arguments ..
      COMPLEX*16         Z
*     ..
*
*  Purpose
*  =======
*
*  ZSLECT returns .TRUE. if the eigenvalue Z is to be selected,
*  otherwise it returns .FALSE.
*  It is used by ZCHK41 to test if ZGEES succesfully sorts eigenvalues,
*  and by ZCHK43 to test if ZGEESX succesfully sorts eigenvalues.
*
*  The common block /SSLCT/ controls how eigenvalues are selected.
*  If SELOPT = 0, then ZSLECT return .TRUE. when real(Z) is less than
*  zero, and .FALSE. otherwise.
*  If SELOPT is at least 1, ZSLECT returns SELVAL(SELOPT) and adds 1
*  to SELOPT, cycling back to 1 at SELMAX.
*
*  Arguments
*  =========
*
*  Z       (input) COMPLEX*16
*          The eigenvalue Z.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   RMIN, X
*     ..
*     .. Scalars in Common ..
      INTEGER            SELDIM, SELOPT
*     ..
*     .. Arrays in Common ..
      LOGICAL            SELVAL( 20 )
      DOUBLE PRECISION   SELWI( 20 ), SELWR( 20 )
*     ..
*     .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX
*     ..
*     .. Executable Statements ..
*
      IF( SELOPT.EQ.0 ) THEN
         ZSLECT = ( DBLE( Z ).LT.ZERO )
      ELSE
         RMIN = ABS( Z-DCMPLX( SELWR( 1 ), SELWI( 1 ) ) )
         ZSLECT = SELVAL( 1 )
         DO 10 I = 2, SELDIM
            X = ABS( Z-DCMPLX( SELWR( I ), SELWI( I ) ) )
            IF( X.LE.RMIN ) THEN
               RMIN = X
               ZSLECT = SELVAL( I )
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of ZSLECT
*
      END
