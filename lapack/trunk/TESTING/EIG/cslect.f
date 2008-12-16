      LOGICAL          FUNCTION CSLECT( Z )
*
*  -- LAPACK test routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     February 2007
*
*     .. Scalar Arguments ..
      COMPLEX            Z
*     ..
*
*  Purpose
*  =======
*
*  CSLECT returns .TRUE. if the eigenvalue Z is to be selected,
*  otherwise it returns .FALSE.
*  It is used by CCHK41 to test if CGEES succesfully sorts eigenvalues,
*  and by CCHK43 to test if CGEESX succesfully sorts eigenvalues.
*
*  The common block /SSLCT/ controls how eigenvalues are selected.
*  If SELOPT = 0, then CSLECT return .TRUE. when real(Z) is less than
*  zero, and .FALSE. otherwise.
*  If SELOPT is at least 1, CSLECT returns SELVAL(SELOPT) and adds 1
*  to SELOPT, cycling back to 1 at SELMAX.
*
*  Arguments
*  =========
*
*  Z       (input) COMPLEX
*          The eigenvalue Z.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               RMIN, X
*     ..
*     .. Scalars in Common ..
      INTEGER            SELDIM, SELOPT
*     ..
*     .. Arrays in Common ..
      LOGICAL            SELVAL( 20 )
      REAL               SELWI( 20 ), SELWR( 20 )
*     ..
*     .. Common blocks ..
      COMMON             / SSLCT / SELOPT, SELDIM, SELVAL, SELWR, SELWI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CMPLX, REAL
*     ..
*     .. Executable Statements ..
*
      IF( SELOPT.EQ.0 ) THEN
         CSLECT = ( REAL( Z ).LT.ZERO )
      ELSE
         RMIN = ABS( Z-CMPLX( SELWR( 1 ), SELWI( 1 ) ) )
         CSLECT = SELVAL( 1 )
         DO 10 I = 2, SELDIM
            X = ABS( Z-CMPLX( SELWR( I ), SELWI( I ) ) )
            IF( X.LE.RMIN ) THEN
               RMIN = X
               CSLECT = SELVAL( I )
            END IF
   10    CONTINUE
      END IF
      RETURN
*
*     End of CSLECT
*
      END
