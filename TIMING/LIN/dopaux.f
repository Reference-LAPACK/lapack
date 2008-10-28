      DOUBLE PRECISION FUNCTION DOPAUX( SUBNAM, M, N, KL, KU, NB )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*(*)       SUBNAM
      INTEGER            KL, KU, M, N, NB
*     ..
*
*  Purpose
*  =======
*
*  DOPAUX computes an approximation of the number of floating point
*  operations used by the subroutine SUBNAM with the given values
*  of the parameters M, N, KL, KU, and NB.
*
*  This version counts operations for the LAPACK auxiliary routines.
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
*          The number of columns of the coefficient matrix.
*          If the matrix is square (such as in a solve routine) then
*          N is the number of right hand sides.  N >= 0.
*
*  KL      (input) INTEGER
*          The lower band width of the coefficient matrix.
*          If needed, 0 <= KL <= M-1.
*
*  KU      (input) INTEGER
*          The upper band width of the coefficient matrix.
*          If needed, 0 <= KU <= N-1.
*
*  NB      (input) INTEGER
*          The block size.  If needed, NB >= 1.
*
*  =====================================================================
*
*     .. Local Scalars ..
      CHARACTER          C1
      CHARACTER*2        C2
      CHARACTER*3        C3
      DOUBLE PRECISION   ADDFAC, ADDS, EK, EM, EN, ENB, MULFAC, MULTS
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      EXTERNAL           LSAME, LSAMEN
*     ..
*     .. Executable Statements ..
*
      DOPAUX = 0
      MULTS = 0
      ADDS = 0
      C1 = SUBNAM( 1: 1 )
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      IF( M.LE.0 .OR. .NOT.( LSAME( C1, 'S' ) .OR. LSAME( C1,
     $    'D' ) .OR. LSAME( C1, 'C' ) .OR. LSAME( C1, 'Z' ) ) ) THEN
         RETURN
      END IF
      IF( LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' ) ) THEN
         MULFAC = 1
         ADDFAC = 1
      ELSE
         MULFAC = 6
         ADDFAC = 2
      END IF
      EM = M
      EN = N
      ENB = NB
*
      IF( LSAMEN( 2, C2, 'LA' ) ) THEN
*
*        xLAULM:  N  =>  M
*
         IF( LSAMEN( 3, C3, 'ULM' ) .OR. LSAMEN( 3, C3, 'UL2' ) ) THEN
            MULTS = ( 1.D0 / 3.D0 )*EM*( -1.D0+EM*EM )
            ADDS = EM*( 1.D0 / 6.D0+EM*( -1.D0 / 2.D0+EM*( 1.D0 /
     $             3.D0 ) ) )
*
*        xLAUUM:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'UUM' ) .OR. LSAMEN( 3, C3, 'UU2' ) )
     $             THEN
            MULTS = EM*( 1.D0 / 3.D0+EM*( 1.D0 / 2.D0+EM*( 1.D0 /
     $              6.D0 ) ) )
            ADDS = ( 1.D0 / 6.D0 )*EM*( -1.D0+EM*EM )
*
*        xLACON:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'CON' ) ) THEN
            MULTS = 3.D0*EM + 3.D0
            ADDS = 4.D0*EM - 3.D0
*
*        xLARF:  M, N  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'RF ' ) ) THEN
            MULTS = 2.D0*EM*EN + EN
            ADDS = 2.D0*EM*EN
*
*        xLARFB:  M, N, SIDE, NB  =>  M, N, KL, NB
*           where KL <= 0 indicates SIDE = 'L'
*           and   KL > 0  indicates SIDE = 'R'
*
         ELSE IF( LSAMEN( 3, C3, 'RFB' ) ) THEN
*
*           KL <= 0:  Code requiring local array
*
            IF( KL.LE.0 ) THEN
               MULTS = EN*ENB*( 2.D0*EM+( ENB+1.D0 ) / 2.D0 )
               ADDS = EN*ENB*( 2.D0*EM+( ENB-1.D0 ) / 2.D0 )
*
*           KL > 0:  Code not requiring local array
*
            ELSE
               MULTS = EN*ENB*( 2.D0*EM+( -ENB / 2.D0+5.D0 / 2.D0 ) )
               ADDS = EN*ENB*( 2.D0*EM+( -ENB / 2.D0-1.D0 / 2.D0 ) )
            END IF
*
*        xLARFG:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'RFG' ) ) THEN
            MULTS = 2.D0*EM + 4.D0
            ADDS = EM + 1.D0
*
*        xLARFT:  M, NB  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'RFT' ) ) THEN
            MULTS = EN*( ( -5.D0 / 6.D0+EN*( 1.D0+EN*( -1.D0 /
     $              6.D0 ) ) )+( EM / 2.D0 )*( EN-1.D0 ) )
            ADDS = EN*( ( 1.D0 / 6.D0 )*( 1.D0-EN*EN )+( EM / 2.D0 )*
     $             ( EN-1.D0 ) )
*
*        xLATRD:  N, K  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'TRD' ) ) THEN
            EK = N
            MULTS = EK*( ( 25.D0 / 6.D0-EK*( 3.D0 / 2.D0+( 5.D0 /
     $              3.D0 )*EK ) )+EM*( 2.D0+2.D0*EK+EM ) )
            ADDS = EK*( ( -1.D0 / 3.D0-( 5.D0 / 3.D0 )*EK*EK )+EM*
     $             ( -1.D0+2.D0*EK+EM ) )
         END IF
*
      END IF
*
      DOPAUX = MULFAC*MULTS + ADDFAC*ADDS
*
      RETURN
*
*     End of DOPAUX
*
      END
