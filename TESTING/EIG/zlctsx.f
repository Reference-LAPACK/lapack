      LOGICAL          FUNCTION ZLCTSX( ALPHA, BETA )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      COMPLEX*16         ALPHA, BETA
*     ..
*
*  Purpose
*  =======
*
*  This function is used to determine what eigenvalues will be
*  selected.  If this is part of the test driver ZDRGSX, do not
*  change the code UNLESS you are testing input examples and not
*  using the built-in examples.
*
*  Arguments
*  =========
*
*  ALPHA   (input) COMPLEX*16
*  BETA    (input) COMPLEX*16
*          parameters to decide whether the pair (ALPHA, BETA) is
*          selected.
*
*  =====================================================================
*
*     .. Parameters ..
*     DOUBLE PRECISION               ZERO
*     PARAMETER          ( ZERO = 0.0E+0 )
*     COMPLEX*16            CZERO
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
            ZLCTSX = .FALSE.
         ELSE
            ZLCTSX = .TRUE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .FALSE.
            I = 0
         END IF
      ELSE
         I = I + 1
         IF( I.LE.N ) THEN
            ZLCTSX = .TRUE.
         ELSE
            ZLCTSX = .FALSE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .TRUE.
            I = 0
         END IF
      END IF
*
*      IF( BETA.EQ.CZERO ) THEN
*         ZLCTSX = ( DBLE( ALPHA ).GT.ZERO )
*      ELSE
*         ZLCTSX = ( DBLE( ALPHA/BETA ).GT.ZERO )
*      END IF
*
      RETURN
*
*     End of ZLCTSX
*
      END
