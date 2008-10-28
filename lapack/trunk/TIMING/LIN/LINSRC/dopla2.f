      DOUBLE PRECISION FUNCTION DOPLA2( SUBNAM, OPTS, M, N, K, L, NB )
*
*  -- LAPACK timing routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*(*)       SUBNAM
      CHARACTER*( * )    OPTS
      INTEGER            K, L, M, N, NB
*     ..
*
*  Purpose
*  =======
*
*  DOPLA2 computes an approximation of the number of floating point
*  operations used by the subroutine SUBNAM with character options
*  OPTS and parameters M, N, K, L, and NB.
*
*  This version counts operations for the LAPACK subroutines that
*  call other LAPACK routines.
*
*  Arguments
*  =========
*
*  SUBNAM  (input) CHARACTER*(*)
*          The name of the subroutine.
*
*  OPTS    (input) CHRACTER*(*)
*          A string of character options to subroutine SUBNAM.
*
*  M       (input) INTEGER
*          The number of rows of the coefficient matrix.
*
*  N       (input) INTEGER
*          The number of columns of the coefficient matrix.
*
*  K       (input) INTEGER
*          A third problem dimension, if needed.
*
*  L       (input) INTEGER
*          A fourth problem dimension, if needed.
*
*  NB      (input) INTEGER
*          The block size.  If needed, NB >= 1.
*
*  Notes
*  =====
*
*  In the comments below, the association is given between arguments
*  in the requested subroutine and local arguments.  For example,
*
*  xORMBR:  VECT // SIDE // TRANS, M, N, K   =>  OPTS, M, N, K
*
*  means that the character string VECT // SIDE // TRANS is passed to
*  the argument OPTS, and the integer parameters M, N, and K are passed
*  to the arguments M, N, and K,
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CORZ, SORD
      CHARACTER          C1, SIDE, UPLO, VECT
      CHARACTER*2        C2
      CHARACTER*3        C3
      CHARACTER(32)      SUB2
      INTEGER            IHI, ILO, ISIDE, MI, NI, NQ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      DOUBLE PRECISION   DOPLA
      EXTERNAL           LSAME, LSAMEN, DOPLA
*     ..
*     .. Executable Statements ..
*
*     ---------------------------------------------------------
*     Initialize DOPLA2 to 0 and do a quick return if possible.
*     ---------------------------------------------------------
*
      DOPLA2 = 0
      C1 = SUBNAM( 1: 1 )
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      SORD = LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' )
      CORZ = LSAME( C1, 'C' ) .OR. LSAME( C1, 'Z' )
      IF( M.LE.0 .OR. .NOT.( SORD .OR. CORZ ) )
     $   RETURN
*
*     -------------------
*     Orthogonal matrices
*     -------------------
*
      IF( ( SORD .AND. LSAMEN( 2, C2, 'OR' ) ) .OR.
     $    ( CORZ .AND. LSAMEN( 2, C2, 'UN' ) ) ) THEN
*
         IF( LSAMEN( 3, C3, 'GBR' ) ) THEN
*
*           -GBR:  VECT, M, N, K  =>  OPTS, M, N, K
*
            VECT = OPTS( 1: 1 )
            IF( LSAME( VECT, 'Q' ) ) THEN
               SUB2 = SUBNAM( 1: 3 ) // 'GQR'
               IF( M.GE.K ) THEN
                  DOPLA2 = DOPLA( SUB2, M, N, K, 0, NB )
               ELSE
                  DOPLA2 = DOPLA( SUB2, M-1, M-1, M-1, 0, NB )
               END IF
            ELSE
               SUB2 = SUBNAM( 1: 3 ) // 'GLQ'
               IF( K.LT.N ) THEN
                  DOPLA2 = DOPLA( SUB2, M, N, K, 0, NB )
               ELSE
                  DOPLA2 = DOPLA( SUB2, N-1, N-1, N-1, 0, NB )
               END IF
            END IF
*
         ELSE IF( LSAMEN( 3, C3, 'MBR' ) ) THEN
*
*           -MBR:  VECT // SIDE // TRANS, M, N, K  =>  OPTS, M, N, K
*
            VECT = OPTS( 1: 1 )
            SIDE = OPTS( 2: 2 )
            IF( LSAME( SIDE, 'L' ) ) THEN
               NQ = M
               ISIDE = 0
            ELSE
               NQ = N
               ISIDE = 1
            END IF
            IF( LSAME( VECT, 'Q' ) ) THEN
               SUB2 = SUBNAM( 1: 3 ) // 'MQR'
               IF( NQ.GE.K ) THEN
                  DOPLA2 = DOPLA( SUB2, M, N, K, ISIDE, NB )
               ELSE IF( ISIDE.EQ.0 ) THEN
                  DOPLA2 = DOPLA( SUB2, M-1, N, NQ-1, ISIDE, NB )
               ELSE
                  DOPLA2 = DOPLA( SUB2, M, N-1, NQ-1, ISIDE, NB )
               END IF
            ELSE
               SUB2 = SUBNAM( 1: 3 ) // 'MLQ'
               IF( NQ.GT.K ) THEN
                  DOPLA2 = DOPLA( SUB2, M, N, K, ISIDE, NB )
               ELSE IF( ISIDE.EQ.0 ) THEN
                  DOPLA2 = DOPLA( SUB2, M-1, N, NQ-1, ISIDE, NB )
               ELSE
                  DOPLA2 = DOPLA( SUB2, M, N-1, NQ-1, ISIDE, NB )
               END IF
            END IF
*
         ELSE IF( LSAMEN( 3, C3, 'GHR' ) ) THEN
*
*           -GHR:  N, ILO, IHI  =>  M, N, K
*
            ILO = N
            IHI = K
            SUB2 = SUBNAM( 1: 3 ) // 'GQR'
            DOPLA2 = DOPLA( SUB2, IHI-ILO, IHI-ILO, IHI-ILO, 0, NB )
*
         ELSE IF( LSAMEN( 3, C3, 'MHR' ) ) THEN
*
*           -MHR:  SIDE // TRANS, M, N, ILO, IHI  =>  OPTS, M, N, K, L
*
            SIDE = OPTS( 1: 1 )
            ILO = K
            IHI = L
            IF( LSAME( SIDE, 'L' ) ) THEN
               MI = IHI - ILO
               NI = N
               ISIDE = -1
            ELSE
               MI = M
               NI = IHI - ILO
               ISIDE = 1
            END IF
            SUB2 = SUBNAM( 1: 3 ) // 'MQR'
            DOPLA2 = DOPLA( SUB2, MI, NI, IHI-ILO, ISIDE, NB )
*
         ELSE IF( LSAMEN( 3, C3, 'GTR' ) ) THEN
*
*           -GTR:  UPLO, N  =>  OPTS, M
*
            UPLO = OPTS( 1: 1 )
            IF( LSAME( UPLO, 'U' ) ) THEN
               SUB2 = SUBNAM( 1: 3 ) // 'GQL'
               DOPLA2 = DOPLA( SUB2, M-1, M-1, M-1, 0, NB )
            ELSE
               SUB2 = SUBNAM( 1: 3 ) // 'GQR'
               DOPLA2 = DOPLA( SUB2, M-1, M-1, M-1, 0, NB )
            END IF
*
         ELSE IF( LSAMEN( 3, C3, 'MTR' ) ) THEN
*
*           -MTR:  SIDE // UPLO // TRANS, M, N  =>  OPTS, M, N
*
            SIDE = OPTS( 1: 1 )
            UPLO = OPTS( 2: 2 )
            IF( LSAME( SIDE, 'L' ) ) THEN
               MI = M - 1
               NI = N
               NQ = M
               ISIDE = -1
            ELSE
               MI = M
               NI = N - 1
               NQ = N
               ISIDE = 1
            END IF
*
            IF( LSAME( UPLO, 'U' ) ) THEN
               SUB2 = SUBNAM( 1: 3 ) // 'MQL'
               DOPLA2 = DOPLA( SUB2, MI, NI, NQ-1, ISIDE, NB )
            ELSE
               SUB2 = SUBNAM( 1: 3 ) // 'MQR'
               DOPLA2 = DOPLA( SUB2, MI, NI, NQ-1, ISIDE, NB )
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of DOPLA2
*
      END
