      SUBROUTINE DLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB,
     $                   ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB,
     $                   DISTA, DISTB )
*
*  -- LAPACK test routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DISTA, DISTB, TYPE
      CHARACTER*3        PATH
      INTEGER            IMAT, KLA, KLB, KUA, KUB, M, MODEA, MODEB, N, P
      DOUBLE PRECISION   ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*
*  Purpose
*  =======
*
*  DLATB9 sets parameters for the matrix generator based on the type of
*  matrix to be generated.
*
*  Arguments
*  =========
*
*  PATH    (input) CHARACTER*3
*          The LAPACK path name.
*
*  IMAT    (input) INTEGER
*          An integer key describing which matrix to generate for this
*          path.
*
*  M       (input) INTEGER
*          The number of rows in the matrix to be generated.
*
*  N       (input) INTEGER
*          The number of columns in the matrix to be generated.
*
*  TYPE    (output) CHARACTER*1
*          The type of the matrix to be generated:
*          = 'S':  symmetric matrix;
*          = 'P':  symmetric positive (semi)definite matrix;
*          = 'N':  nonsymmetric matrix.
*
*  KL      (output) INTEGER
*          The lower band width of the matrix to be generated.
*
*  KU      (output) INTEGER
*          The upper band width of the matrix to be generated.
*
*  ANORM   (output) DOUBLE PRECISION
*          The desired norm of the matrix to be generated.  The diagonal
*          matrix of singular values or eigenvalues is scaled by this
*          value.
*
*  MODE    (output) INTEGER
*          A key indicating how to choose the vector of eigenvalues.
*
*  CNDNUM  (output) DOUBLE PRECISION
*          The desired condition number.
*
*  DIST    (output) CHARACTER*1
*          The type of distribution to be used by the random number
*          generator.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   SHRINK, TENTH
      PARAMETER          ( SHRINK = 0.25D0, TENTH = 0.1D+0 )
      DOUBLE PRECISION   ONE, TEN
      PARAMETER          ( ONE = 1.0D+0, TEN = 1.0D+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST
      DOUBLE PRECISION   BADC1, BADC2, EPS, LARGE, SMALL
*     ..
*     .. External Functions ..
      LOGICAL            LSAMEN
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAMEN, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD
*     ..
*     .. Save statement ..
      SAVE               EPS, SMALL, LARGE, BADC1, BADC2, FIRST
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
*     Set some constants for use in the subroutine.
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         EPS = DLAMCH( 'Precision' )
         BADC2 = TENTH / EPS
         BADC1 = SQRT( BADC2 )
         SMALL = DLAMCH( 'Safe minimum' )
         LARGE = ONE / SMALL
*
*        If it looks like we're on a Cray, take the square root of
*        SMALL and LARGE to avoid overflow and underflow problems.
*
         CALL DLABAD( SMALL, LARGE )
         SMALL = SHRINK*( SMALL / EPS )
         LARGE = ONE / SMALL
      END IF
*
*     Set some parameters we don't plan to change.
*
      TYPE = 'N'
      DISTA = 'S'
      DISTB = 'S'
      MODEA = 3
      MODEB = 4
*
*     Set the lower and upper bandwidths.
*
      IF( LSAMEN( 3, PATH, 'GRQ' ) .OR. LSAMEN( 3, PATH, 'LSE' ) .OR.
     $    LSAMEN( 3, PATH, 'GSV' ) ) THEN
*
*        A: M by N, B: P by N
*
         IF( IMAT.EQ.1 ) THEN
*
*           A: diagonal, B: upper triangular
*
            KLA = 0
            KUA = 0
            KLB = 0
            KUB = MAX( N-1, 0 )
*
         ELSE IF( IMAT.EQ.2 ) THEN
*
*           A: upper triangular, B: upper triangular
*
            KLA = 0
            KUA = MAX( N-1, 0 )
            KLB = 0
            KUB = MAX( N-1, 0 )
*
         ELSE IF( IMAT.EQ.3 ) THEN
*
*           A: lower triangular, B: upper triangular
*
            KLA = MAX( M-1, 0 )
            KUA = 0
            KLB = 0
            KUB = MAX( N-1, 0 )
*
         ELSE
*
*           A: general dense, B: general dense
*
            KLA = MAX( M-1, 0 )
            KUA = MAX( N-1, 0 )
            KLB = MAX( P-1, 0 )
            KUB = MAX( N-1, 0 )
*
         END IF
*
      ELSE IF( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GLM' ) )
     $          THEN
*
*        A: N by M, B: N by P
*
         IF( IMAT.EQ.1 ) THEN
*
*           A: diagonal, B: lower triangular
*
            KLA = 0
            KUA = 0
            KLB = MAX( N-1, 0 )
            KUB = 0
         ELSE IF( IMAT.EQ.2 ) THEN
*
*           A: lower triangular, B: diagonal
*
            KLA = MAX( N-1, 0 )
            KUA = 0
            KLB = 0
            KUB = 0
*
         ELSE IF( IMAT.EQ.3 ) THEN
*
*           A: lower triangular, B: upper triangular
*
            KLA = MAX( N-1, 0 )
            KUA = 0
            KLB = 0
            KUB = MAX( P-1, 0 )
*
         ELSE
*
*           A: general dense, B: general dense
*
            KLA = MAX( N-1, 0 )
            KUA = MAX( M-1, 0 )
            KLB = MAX( N-1, 0 )
            KUB = MAX( P-1, 0 )
         END IF
*
      END IF
*
*     Set the condition number and norm.
*
      CNDNMA = TEN*TEN
      CNDNMB = TEN
      IF( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' ) .OR.
     $    LSAMEN( 3, PATH, 'GSV' ) ) THEN
         IF( IMAT.EQ.5 ) THEN
            CNDNMA = BADC1
            CNDNMB = BADC1
         ELSE IF( IMAT.EQ.6 ) THEN
            CNDNMA = BADC2
            CNDNMB = BADC2
         ELSE IF( IMAT.EQ.7 ) THEN
            CNDNMA = BADC1
            CNDNMB = BADC2
         ELSE IF( IMAT.EQ.8 ) THEN
            CNDNMA = BADC2
            CNDNMB = BADC1
         END IF
      END IF
*
      ANORM = TEN
      BNORM = TEN*TEN*TEN
      IF( LSAMEN( 3, PATH, 'GQR' ) .OR. LSAMEN( 3, PATH, 'GRQ' ) ) THEN
         IF( IMAT.EQ.7 ) THEN
            ANORM = SMALL
            BNORM = LARGE
         ELSE IF( IMAT.EQ.8 ) THEN
            ANORM = LARGE
            BNORM = SMALL
         END IF
      END IF
*
      IF( N.LE.1 ) THEN
         CNDNMA = ONE
         CNDNMB = ONE
      END IF
*
      RETURN
*
*     End of DLATB9
*
      END
