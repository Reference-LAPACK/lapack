      DOUBLE PRECISION FUNCTION DOPGB( SUBNAM, M, N, KL, KU, IPIV )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*(*)       SUBNAM
      INTEGER            KL, KU, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
*     ..
*
*  Purpose
*  =======
*
*  DOPGB counts operations for the LU factorization of a band matrix
*  xGBTRF.
*
*  Arguments
*  =========
*
*  SUBNAM  (input) CHARACTER*(*)
*          The name of the subroutine.
*
*  M       (input) INTEGER
*          The number of rows of the coefficient matrix.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the coefficient matrix.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals of the matrix.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals of the matrix.  KU >= 0.
*
*  IPIV    (input)  INTEGER array, dimension (min(M,N))
*          The vector of pivot indices from DGBTRF or ZGBTRF.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CORZ, SORD
      CHARACTER          C1
      CHARACTER*2        C2
      CHARACTER*3        C3
      INTEGER            I, J, JP, JU, KM
      DOUBLE PRECISION   ADDFAC, ADDS, MULFAC, MULTS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           LSAME, LSAMEN
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
      DOPGB = 0
      MULTS = 0
      ADDS = 0
      C1 = SUBNAM( 1: 1 )
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      SORD = LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' )
      CORZ = LSAME( C1, 'C' ) .OR. LSAME( C1, 'Z' )
      IF( .NOT.( SORD .OR. CORZ ) )
     $   RETURN
      IF( LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' ) ) THEN
         ADDFAC = 1
         MULFAC = 1
      ELSE
         ADDFAC = 2
         MULFAC = 6
      END IF
*
*     --------------------------
*     GB:  General Band matrices
*     --------------------------
*
      IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
*        xGBTRF:  M, N, KL, KU  =>  M, N, KL, KU
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            JU = 1
            DO 10 J = 1, MIN( M, N )
               KM = MIN( KL, M-J )
               JP = IPIV( J )
               JU = MAX( JU, MIN( JP+KU, N ) )
               IF( KM.GT.0 ) THEN
                  MULTS = MULTS + KM*( 1+JU-J )
                  ADDS = ADDS + KM*( JU-J )
               END IF
   10       CONTINUE
         END IF
*
*     ---------------------------------
*     GT:  General Tridiagonal matrices
*     ---------------------------------
*
      ELSE IF( LSAMEN( 2, C2, 'GT' ) ) THEN
*
*        xGTTRF:  N  =>  M
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            MULTS = 2*( M-1 )
            ADDS = M - 1
            DO 20 I = 1, M - 2
               IF( IPIV( I ).NE.I )
     $            MULTS = MULTS + 1
   20       CONTINUE
*
*        xGTTRS:  N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = 4*N*( M-1 )
            ADDS = 3*N*( M-1 )
*
*        xGTSV:   N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'SV ' ) ) THEN
            MULTS = ( 4*N+2 )*( M-1 )
            ADDS = ( 3*N+1 )*( M-1 )
            DO 30 I = 1, M - 2
               IF( IPIV( I ).NE.I )
     $            MULTS = MULTS + 1
   30       CONTINUE
         END IF
      END IF
*
      DOPGB = MULFAC*MULTS + ADDFAC*ADDS
      RETURN
*
*     End of DOPGB
*
      END
