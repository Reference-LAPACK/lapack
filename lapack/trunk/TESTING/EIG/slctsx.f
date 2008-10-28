      LOGICAL          FUNCTION SLCTSX( AR, AI, BETA )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      REAL               AI, AR, BETA
*     ..
*
*  Purpose
*  =======
*
*  This function is used to determine what eigenvalues will be
*  selected.  If this is part of the test driver SDRGSX, do not
*  change the code UNLESS you are testing input examples and not
*  using the built-in examples.
*
*  Arguments
*  =========
*
*  AR      (input) REAL
*          The numerator of the real part of a complex eigenvalue
*          (AR/BETA) + i*(AI/BETA).
*
*  AI      (input) REAL
*          The numerator of the imaginary part of a complex eigenvalue
*          (AR/BETA) + i*(AI).
*
*  BETA    (input) REAL
*          The denominator part of a complex eigenvalue
*          (AR/BETA) + i*(AI/BETA).
*
*  =====================================================================
*
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
            SLCTSX = .FALSE.
         ELSE
            SLCTSX = .TRUE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .FALSE.
            I = 0
         END IF
      ELSE
         I = I + 1
         IF( I.LE.N ) THEN
            SLCTSX = .TRUE.
         ELSE
            SLCTSX = .FALSE.
         END IF
         IF( I.EQ.MPLUSN ) THEN
            FS = .TRUE.
            I = 0
         END IF
      END IF
*
*       IF( AR/BETA.GT.0.0 )THEN
*          SLCTSX = .TRUE.
*       ELSE
*          SLCTSX = .FALSE.
*       END IF
*
      RETURN
*
*     End of SLCTSX
*
      END
