      LOGICAL          FUNCTION CLCTSX( ALPHA, BETA )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      COMPLEX            ALPHA, BETA
*     ..
*
*  Purpose
*  =======
*
*  This function is used to determine what eigenvalues will be
*  selected.  If this is part of the test driver CDRGSX, do not
*  change the code UNLESS you are testing input examples and not
*  using the built-in examples.
*
*  Arguments
*  =========
*
*  ALPHA   (input) COMPLEX
*  BETA    (input) COMPLEX
*          parameters to decide whether the pair (ALPHA, BETA) is
*          selected.
*
*  =====================================================================
*
*     .. Parameters ..
*     REAL               ZERO
*     PARAMETER          ( ZERO = 0.0E+0 )
*     COMPLEX            CZERO
*     PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Scalars in Common ..
      LOGICAL            FS
      INTEGER            I, M, MPLUSN, N
*     ..
*     .. Common blocks ..
      COMMON             / MN / M, N, MPLUSN, I, FS
*     ..
*     .. Save statement ..
      SAVE
*     ..
*     .. Executable Statements ..
*
      IF( FS ) THEN
         I = I + 1
         IF( I.LE.M ) THEN
            CLCTSX = .FALSE.
         ELSE
            CLCTSX = .TRUE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .FALSE.
            I = 0
         END IF
      ELSE
         I = I + 1
         IF( I.LE.N ) THEN
            CLCTSX = .TRUE.
         ELSE
            CLCTSX = .FALSE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .TRUE.
            I = 0
         END IF
      END IF
*
*      IF( BETA.EQ.CZERO ) THEN
*         CLCTSX = ( REAL( ALPHA ).GT.ZERO )
*      ELSE
*         CLCTSX = ( REAL( ALPHA/BETA ).GT.ZERO )
*      END IF
*
      RETURN
*
*     End of CLCTSX
*
      END
