      DOUBLE PRECISION FUNCTION DOPLA( SUBNAM, M, N, KL, KU, NB )
*
*  -- LAPACK timing routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER*(*)       SUBNAM
      INTEGER            KL, KU, M, N, NB
*     ..
*
*  Purpose
*  =======
*
*  DOPLA computes an approximation of the number of floating point
*  operations used by the subroutine SUBNAM with the given values
*  of the parameters M, N, KL, KU, and NB.
*
*  This version counts operations for the LAPACK subroutines.
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
*          For solve routine when the matrix is square,
*          N is the number of right hand sides.  N >= 0.
*
*  KL      (input) INTEGER
*          The lower band width of the coefficient matrix.
*          If needed, 0 <= KL <= M-1.
*          For xGEQRS, KL is the number of right hand sides.
*
*  KU      (input) INTEGER
*          The upper band width of the coefficient matrix.
*          If needed, 0 <= KU <= N-1.
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
*  xGETRS:  N, NRHS  =>  M, N
*
*  means that arguments N and NRHS in DGETRS are passed to arguments
*  M and N in this procedure.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CORZ, SORD
      CHARACTER          C1
      CHARACTER*2        C2
      CHARACTER*3        C3
      INTEGER            I
      DOUBLE PRECISION   ADDFAC, ADDS, EK, EM, EMN, EN, MULFAC, MULTS,
     $                   WL, WU
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
*     --------------------------------------------------------
*     Initialize DOPLA to 0 and do a quick return if possible.
*     --------------------------------------------------------
*
      DOPLA = 0
      MULTS = 0
      ADDS = 0
      C1 = SUBNAM( 1: 1 )
      C2 = SUBNAM( 2: 3 )
      C3 = SUBNAM( 4: 6 )
      SORD = LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' )
      CORZ = LSAME( C1, 'C' ) .OR. LSAME( C1, 'Z' )
      IF( M.LE.0 .OR. .NOT.( SORD .OR. CORZ ) )
     $   RETURN
*
*     ---------------------------------------------------------
*     If the coefficient matrix is real, count each add as 1
*     operation and each multiply as 1 operation.
*     If the coefficient matrix is complex, count each add as 2
*     operations and each multiply as 6 operations.
*     ---------------------------------------------------------
*
      IF( LSAME( C1, 'S' ) .OR. LSAME( C1, 'D' ) ) THEN
         ADDFAC = 1
         MULFAC = 1
      ELSE
         ADDFAC = 2
         MULFAC = 6
      END IF
      EM = M
      EN = N
      EK = KL
*
*     ---------------------------------
*     GE:  GEneral rectangular matrices
*     ---------------------------------
*
      IF( LSAMEN( 2, C2, 'GE' ) ) THEN
*
*        xGETRF:  M, N  =>  M, N
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            EMN = MIN( M, N )
            ADDS = EMN*( EM*EN-( EM+EN )*( EMN+1.D0 ) / 2.D0+
     $             ( EMN+1.D0 )*( 2.D0*EMN+1.D0 ) / 6.D0 )
            MULTS = ADDS + EMN*( EM-( EMN+1.D0 ) / 2.D0 )
*
*        xGETRS:  N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*EM*EM
            ADDS = EN*( EM*( EM-1.D0 ) )
*
*        xGETRI:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            MULTS = EM*( 5.D0 / 6.D0+EM*( 1.D0 / 2.D0+EM*( 2.D0 /
     $              3.D0 ) ) )
            ADDS = EM*( 5.D0 / 6.D0+EM*( -3.D0 / 2.D0+EM*( 2.D0 /
     $             3.D0 ) ) )
*
*        xGEQRF or xGEQLF:  M, N  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'QRF' ) .OR.
     $            LSAMEN( 3, C3, 'QR2' ) .OR.
     $            LSAMEN( 3, C3, 'QLF' ) .OR. LSAMEN( 3, C3, 'QL2' ) )
     $             THEN
            IF( M.GE.N ) THEN
               MULTS = EN*( ( ( 23.D0 / 6.D0 )+EM+EN / 2.D0 )+EN*
     $                 ( EM-EN / 3.D0 ) )
               ADDS = EN*( ( 5.D0 / 6.D0 )+EN*
     $                ( 1.D0 / 2.D0+( EM-EN / 3.D0 ) ) )
            ELSE
               MULTS = EM*( ( ( 23.D0 / 6.D0 )+2.D0*EN-EM / 2.D0 )+EM*
     $                 ( EN-EM / 3.D0 ) )
               ADDS = EM*( ( 5.D0 / 6.D0 )+EN-EM / 2.D0+EM*
     $                ( EN-EM / 3.D0 ) )
            END IF
*
*        xGERQF or xGELQF:  M, N  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'RQF' ) .OR.
     $            LSAMEN( 3, C3, 'RQ2' ) .OR.
     $            LSAMEN( 3, C3, 'LQF' ) .OR. LSAMEN( 3, C3, 'LQ2' ) )
     $             THEN
            IF( M.GE.N ) THEN
               MULTS = EN*( ( ( 29.D0 / 6.D0 )+EM+EN / 2.D0 )+EN*
     $                 ( EM-EN / 3.D0 ) )
               ADDS = EN*( ( 5.D0 / 6.D0 )+EM+EN*
     $                ( -1.D0 / 2.D0+( EM-EN / 3.D0 ) ) )
            ELSE
               MULTS = EM*( ( ( 29.D0 / 6.D0 )+2.D0*EN-EM / 2.D0 )+EM*
     $                 ( EN-EM / 3.D0 ) )
               ADDS = EM*( ( 5.D0 / 6.D0 )+EM / 2.D0+EM*
     $                ( EN-EM / 3.D0 ) )
            END IF
*
*        xGEQPF: M, N => M, N
*
         ELSE IF( LSAMEN( 3, C3, 'QPF' ) ) THEN
            EMN = MIN( M, N )
            MULTS = 2*EN*EN + EMN*( 3*EM+5*EN+2*EM*EN-( EMN+1 )*
     $              ( 4+EN+EM-( 2*EMN+1 ) / 3 ) )
            ADDS = EN*EN + EMN*( 2*EM+EN+2*EM*EN-( EMN+1 )*
     $             ( 2+EN+EM-( 2*EMN+1 ) / 3 ) )
*
*        xGEQRS or xGERQS:  M, N, NRHS  =>  M, N, KL
*
         ELSE IF( LSAMEN( 3, C3, 'QRS' ) .OR. LSAMEN( 3, C3, 'RQS' ) )
     $             THEN
            MULTS = EK*( EN*( 2.D0-EK )+EM*
     $              ( 2.D0*EN+( EM+1.D0 ) / 2.D0 ) )
            ADDS = EK*( EN*( 1.D0-EK )+EM*
     $             ( 2.D0*EN+( EM-1.D0 ) / 2.D0 ) )
*
*        xGELQS or xGEQLS:  M, N, NRHS  =>  M, N, KL
*
         ELSE IF( LSAMEN( 3, C3, 'LQS' ) .OR. LSAMEN( 3, C3, 'QLS' ) )
     $             THEN
            MULTS = EK*( EM*( 2.D0-EK )+EN*
     $              ( 2.D0*EM+( EN+1.D0 ) / 2.D0 ) )
            ADDS = EK*( EM*( 1.D0-EK )+EN*
     $             ( 2.D0*EM+( EN-1.D0 ) / 2.D0 ) )
*
*        xGEBRD:  M, N  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'BRD' ) ) THEN
            IF( M.GE.N ) THEN
               MULTS = EN*( 20.D0 / 3.D0+EN*
     $                 ( 2.D0+( 2.D0*EM-( 2.D0 / 3.D0 )*EN ) ) )
               ADDS = EN*( 5.D0 / 3.D0+( EN-EM )+EN*
     $                ( 2.D0*EM-( 2.D0 / 3.D0 )*EN ) )
            ELSE
               MULTS = EM*( 20.D0 / 3.D0+EM*
     $                 ( 2.D0+( 2.D0*EN-( 2.D0 / 3.D0 )*EM ) ) )
               ADDS = EM*( 5.D0 / 3.D0+( EM-EN )+EM*
     $                ( 2.D0*EN-( 2.D0 / 3.D0 )*EM ) )
            END IF
*
*        xGEHRD:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'HRD' ) ) THEN
            IF( M.EQ.1 ) THEN
               MULTS = 0.D0
               ADDS = 0.D0
            ELSE
               MULTS = -13.D0 + EM*( -7.D0 / 6.D0+EM*
     $                 ( 0.5D0+EM*( 5.D0 / 3.D0 ) ) )
               ADDS = -8.D0 + EM*( -2.D0 / 3.D0+EM*
     $                ( -1.D0+EM*( 5.D0 / 3.D0 ) ) )
            END IF
*
         END IF
*
*     ----------------------------
*     GB:  General Banded matrices
*     ----------------------------
*        Note:  The operation count is overestimated because
*        it is assumed that the factor U fills in to the maximum
*        extent, i.e., that its bandwidth goes from KU to KL + KU.
*
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
*        xGBTRF:  M, N, KL, KU  =>  M, N, KL, KU
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            DO 10 I = MIN( M, N ), 1, -1
               WL = MAX( 0, MIN( KL, M-I ) )
               WU = MAX( 0, MIN( KL+KU, N-I ) )
               MULTS = MULTS + WL*( 1.D0+WU )
               ADDS = ADDS + WL*WU
   10       CONTINUE
*
*        xGBTRS:  N, NRHS, KL, KU  =>  M, N, KL, KU
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            WL = MAX( 0, MIN( KL, M-1 ) )
            WU = MAX( 0, MIN( KL+KU, M-1 ) )
            MULTS = EN*( EM*( WL+1.D0+WU )-0.5D0*
     $              ( WL*( WL+1.D0 )+WU*( WU+1.D0 ) ) )
            ADDS = EN*( EM*( WL+WU )-0.5D0*
     $             ( WL*( WL+1.D0 )+WU*( WU+1.D0 ) ) )
*
         END IF
*
*     --------------------------------------
*     PO:  POsitive definite matrices
*     PP:  Positive definite Packed matrices
*     --------------------------------------
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) .OR. LSAMEN( 2, C2, 'PP' ) ) THEN
*
*        xPOTRF:  N  =>  M
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            MULTS = EM*( 1.D0 / 3.D0+EM*( 1.D0 / 2.D0+EM*( 1.D0 /
     $              6.D0 ) ) )
            ADDS = ( 1.D0 / 6.D0 )*EM*( -1.D0+EM*EM )
*
*        xPOTRS:  N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*( EM*( EM+1.D0 ) )
            ADDS = EN*( EM*( EM-1.D0 ) )
*
*        xPOTRI:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            MULTS = EM*( 2.D0 / 3.D0+EM*( 1.D0+EM*( 1.D0 / 3.D0 ) ) )
            ADDS = EM*( 1.D0 / 6.D0+EM*( -1.D0 / 2.D0+EM*( 1.D0 /
     $             3.D0 ) ) )
*
         END IF
*
*     ------------------------------------
*     PB:  Positive definite Band matrices
*     ------------------------------------
*
      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
*
*        xPBTRF:  N, K  =>  M, KL
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            MULTS = EK*( -2.D0 / 3.D0+EK*( -1.D0+EK*( -1.D0 / 3.D0 ) ) )
     $               + EM*( 1.D0+EK*( 3.D0 / 2.D0+EK*( 1.D0 / 2.D0 ) ) )
            ADDS = EK*( -1.D0 / 6.D0+EK*( -1.D0 / 2.D0+EK*( -1.D0 /
     $             3.D0 ) ) ) + EM*( EK / 2.D0*( 1.D0+EK ) )
*
*        xPBTRS:  N, NRHS, K  =>  M, N, KL
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*( ( 2*EM-EK )*( EK+1.D0 ) )
            ADDS = EN*( EK*( 2*EM-( EK+1.D0 ) ) )
*
         END IF
*
*     ----------------------------------
*     PT:  Positive definite Tridiagonal
*     ----------------------------------
*
      ELSE IF( LSAMEN( 2, C2, 'PT' ) ) THEN
*
*        xPTTRF:  N  =>  M
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            MULTS = 2*( EM-1 )
            ADDS = EM - 1
*
*        xPTTRS:  N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*( 3*EM-2 )
            ADDS = EN*( 2*( EM-1 ) )
*
*        xPTSV:  N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'SV ' ) ) THEN
            MULTS = 2*( EM-1 ) + EN*( 3*EM-2 )
            ADDS = EM - 1 + EN*( 2*( EM-1 ) )
         END IF
*
*     --------------------------------------------------------
*     SY:  SYmmetric indefinite matrices
*     SP:  Symmetric indefinite Packed matrices
*     HE:  HErmitian indefinite matrices (complex only)
*     HP:  Hermitian indefinite Packed matrices (complex only)
*     --------------------------------------------------------
*
      ELSE IF( LSAMEN( 2, C2, 'SY' ) .OR. LSAMEN( 2, C2, 'SP' ) .OR.
     $         LSAMEN( 3, SUBNAM, 'CHE' ) .OR.
     $         LSAMEN( 3, SUBNAM, 'ZHE' ) .OR.
     $         LSAMEN( 3, SUBNAM, 'CHP' ) .OR.
     $         LSAMEN( 3, SUBNAM, 'ZHP' ) ) THEN
*
*        xSYTRF:  N  =>  M
*
         IF( LSAMEN( 3, C3, 'TRF' ) ) THEN
            MULTS = EM*( 10.D0 / 3.D0+EM*
     $              ( 1.D0 / 2.D0+EM*( 1.D0 / 6.D0 ) ) )
            ADDS = EM / 6.D0*( -1.D0+EM*EM )
*
*        xSYTRS:  N, NRHS  =>  M, N
*
         ELSE IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*EM*EM
            ADDS = EN*( EM*( EM-1.D0 ) )
*
*        xSYTRI:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            MULTS = EM*( 2.D0 / 3.D0+EM*EM*( 1.D0 / 3.D0 ) )
            ADDS = EM*( -1.D0 / 3.D0+EM*EM*( 1.D0 / 3.D0 ) )
*
*        xSYTRD, xSYTD2:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'TRD' ) .OR. LSAMEN( 3, C3, 'TD2' ) )
     $             THEN
            IF( M.EQ.1 ) THEN
               MULTS = 0.D0
               ADDS = 0.D0
            ELSE
               MULTS = -15.D0 + EM*( -1.D0 / 6.D0+EM*
     $                 ( 5.D0 / 2.D0+EM*( 2.D0 / 3.D0 ) ) )
               ADDS = -4.D0 + EM*( -8.D0 / 3.D0+EM*
     $                ( 1.D0+EM*( 2.D0 / 3.D0 ) ) )
            END IF
         END IF
*
*     -------------------
*     Triangular matrices
*     -------------------
*
      ELSE IF( LSAMEN( 2, C2, 'TR' ) .OR. LSAMEN( 2, C2, 'TP' ) ) THEN
*
*        xTRTRS:  N, NRHS  =>  M, N
*
         IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*EM*( EM+1.D0 ) / 2.D0
            ADDS = EN*EM*( EM-1.D0 ) / 2.D0
*
*        xTRTRI:  N  =>  M
*
         ELSE IF( LSAMEN( 3, C3, 'TRI' ) ) THEN
            MULTS = EM*( 1.D0 / 3.D0+EM*( 1.D0 / 2.D0+EM*( 1.D0 /
     $              6.D0 ) ) )
            ADDS = EM*( 1.D0 / 3.D0+EM*( -1.D0 / 2.D0+EM*( 1.D0 /
     $             6.D0 ) ) )
*
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
*
*        xTBTRS:  N, NRHS, K  =>  M, N, KL
*
         IF( LSAMEN( 3, C3, 'TRS' ) ) THEN
            MULTS = EN*( EM*( EM+1.D0 ) / 2.D0-( EM-EK-1.D0 )*
     $              ( EM-EK ) / 2.D0 )
            ADDS = EN*( EM*( EM-1.D0 ) / 2.D0-( EM-EK-1.D0 )*( EM-EK ) /
     $             2.D0 )
         END IF
*
*     --------------------
*     Trapezoidal matrices
*     --------------------
*
      ELSE IF( LSAMEN( 2, C2, 'TZ' ) ) THEN
*
*        xTZRQF:  M, N => M, N
*
         IF( LSAMEN( 3, C3, 'RQF' ) ) THEN
            EMN = MIN( M, N )
            MULTS = 3*EM*( EN-EM+1 ) + ( 2*EN-2*EM+3 )*
     $              ( EM*EM-EMN*( EMN+1 ) / 2 )
            ADDS = ( EN-EM+1 )*( EM+2*EM*EM-EMN*( EMN+1 ) )
         END IF
*
*     -------------------
*     Orthogonal matrices
*     -------------------
*
      ELSE IF( ( SORD .AND. LSAMEN( 2, C2, 'OR' ) ) .OR.
     $         ( CORZ .AND. LSAMEN( 2, C2, 'UN' ) ) ) THEN
*
*        -MQR, -MLQ, -MQL, or -MRQ:  M, N, K, SIDE  =>  M, N, KL, KU
*           where KU<= 0 indicates SIDE = 'L'
*           and   KU> 0  indicates SIDE = 'R'
*
         IF( LSAMEN( 3, C3, 'MQR' ) .OR. LSAMEN( 3, C3, 'MLQ' ) .OR.
     $       LSAMEN( 3, C3, 'MQL' ) .OR. LSAMEN( 3, C3, 'MRQ' ) ) THEN
            IF( KU.LE.0 ) THEN
               MULTS = EK*EN*( 2.D0*EM+2.D0-EK )
               ADDS = EK*EN*( 2.D0*EM+1.D0-EK )
            ELSE
               MULTS = EK*( EM*( 2.D0*EN-EK )+
     $                 ( EM+EN+( 1.D0-EK ) / 2.D0 ) )
               ADDS = EK*EM*( 2.D0*EN+1.D0-EK )
            END IF
*
*        -GQR or -GQL:  M, N, K  =>  M, N, KL
*
         ELSE IF( LSAMEN( 3, C3, 'GQR' ) .OR. LSAMEN( 3, C3, 'GQL' ) )
     $             THEN
            MULTS = EK*( -5.D0 / 3.D0+( 2.D0*EN-EK )+
     $              ( 2.D0*EM*EN+EK*( ( 2.D0 / 3.D0 )*EK-EM-EN ) ) )
            ADDS = EK*( 1.D0 / 3.D0+( EN-EM )+
     $             ( 2.D0*EM*EN+EK*( ( 2.D0 / 3.D0 )*EK-EM-EN ) ) )
*
*        -GLQ or -GRQ:  M, N, K  =>  M, N, KL
*
         ELSE IF( LSAMEN( 3, C3, 'GLQ' ) .OR. LSAMEN( 3, C3, 'GRQ' ) )
     $             THEN
            MULTS = EK*( -2.D0 / 3.D0+( EM+EN-EK )+
     $              ( 2.D0*EM*EN+EK*( ( 2.D0 / 3.D0 )*EK-EM-EN ) ) )
            ADDS = EK*( 1.D0 / 3.D0+( EM-EN )+
     $             ( 2.D0*EM*EN+EK*( ( 2.D0 / 3.D0 )*EK-EM-EN ) ) )
*
         END IF
*
      END IF
*
      DOPLA = MULFAC*MULTS + ADDFAC*ADDS
*
      RETURN
*
*     End of DOPLA
*
      END
