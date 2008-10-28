      SUBROUTINE CTIMGT( LINE, NM, MVAL, NNS, NSVAL, NLDA, LDAVAL,
     $                   TIMMIN, A, B, IWORK, RESLTS, LDR1, LDR2, LDR3,
     $                   NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NLDA, NM, NNS, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NSVAL( * )
      REAL               RESLTS( LDR1, LDR2, LDR3, * )
      COMPLEX            A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  CTIMGT times CGTTRF, -TRS, -SV, and -SL.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          The input line that requested this routine.  The first six
*          characters contain either the name of a subroutine or a
*          generic path name.  The remaining characters may be used to
*          specify the individual routines to be timed.  See ATIMIN for
*          a full description of the format of the input line.
*
*  NM      (input) INTEGER
*          The number of values of M contained in the vector MVAL.
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix size M.
*
*  NNS     (input) INTEGER
*          The number of values of NRHS contained in the vector NSVAL.
*
*  NSVAL   (input) INTEGER array, dimension (NNS)
*          The values of the number of right hand sides NRHS.
*
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  TIMMIN  (input) REAL
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) COMPLEX array, dimension (NMAX*4)
*          where NMAX is the maximum value permitted for N.
*
*  B       (workspace) COMPLEX array, dimension (LDAMAX*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (output) REAL array, dimension
*                   (LDR1,LDR2,LDR3,NSUBS+1)
*          The timing results for each subroutine over the relevant
*          values of N.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= 1.
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
*
*  LDR3    (input) INTEGER
*          The third dimension of RESLTS.  LDR3 >= max(1,NLDA).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 4 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANS
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IM, INFO, ISUB, ITRANS, LDB,
     $                   M, N, NRHS
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          TRANSS( 3 )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            LAVAL( 1 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SMFLOP, SOPGB
      EXTERNAL           SECOND, SMFLOP, SOPGB
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, CGTSL, CGTSV, CGTTRF, CGTTRS,
     $                   CTIMMG, SPRTBL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CGTTRF', 'CGTTRS', 'CGTSV ',
     $                   'CGTSL ' /
      DATA               TRANSS / 'N', 'T', 'C' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'GT'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 180
*
*     Check that N <= LDA for the input values.
*
      DO 10 ISUB = 2, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 10
         CNAME = SUBNAM( ISUB )
         CALL ATIMCK( 2, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
         IF( INFO.GT.0 ) THEN
            WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
            TIMSUB( ISUB ) = .FALSE.
         END IF
   10 CONTINUE
*
*     Do for each value of M:
*
      DO 150 IM = 1, NM
*
         M = MVAL( IM )
         N = MAX( M, 1 )
*
*        Time CGTTRF
*
         IF( TIMSUB( 1 ) ) THEN
            CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
            IC = 0
            S1 = SECOND( )
   20       CONTINUE
            CALL CGTTRF( M, A, A( N ), A( 2*N ), A( 3*N-2 ), IWORK,
     $                   INFO )
            S2 = SECOND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
               GO TO 20
            END IF
*
*           Subtract the time used in CTIMMG.
*
            ICL = 1
            S1 = SECOND( )
   30       CONTINUE
            S2 = SECOND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
               GO TO 30
            END IF
*
            TIME = ( TIME-UNTIME ) / REAL( IC )
            OPS = SOPGB( 'CGTTRF', M, M, 1, 1, IWORK )
            RESLTS( 1, IM, 1, 1 ) = SMFLOP( OPS, TIME, INFO )
*
         ELSE IF( TIMSUB( 2 ) ) THEN
            CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
         END IF
*
*        Generate another matrix and factor it using CGTTRF so
*        that the factored form can be used in timing the other
*        routines.
*
         IF( IC.NE.1 )
     $      CALL CGTTRF( M, A, A( N ), A( 2*N ), A( 3*N-2 ), IWORK,
     $                   INFO )
*
*        Time CGTTRS
*
         IF( TIMSUB( 2 ) ) THEN
            DO 80 ITRANS = 1, 3
               TRANS = TRANSS( ITRANS )
               IF( ITRANS.EQ.1 ) THEN
                  ISUB = 2
               ELSE
                  ISUB = ITRANS + 3
               END IF
               DO 70 ILDA = 1, NLDA
                  LDB = LDAVAL( ILDA )
                  DO 60 I = 1, NNS
                     NRHS = NSVAL( I )
                     CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
   40                CONTINUE
                     CALL CGTTRS( TRANS, M, NRHS, A, A( N ), A( 2*N ),
     $                            A( 3*N-2 ), IWORK, B, LDB, INFO )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                        GO TO 40
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
   50                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                        GO TO 50
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPGB( 'CGTTRS', M, NRHS, 0, 0, IWORK )
                     RESLTS( I, IM, ILDA, ISUB ) = SMFLOP( OPS, TIME,
     $                  INFO )
   60             CONTINUE
   70          CONTINUE
   80       CONTINUE
         END IF
*
         IF( TIMSUB( 3 ) ) THEN
            DO 120 ILDA = 1, NLDA
               LDB = LDAVAL( ILDA )
               DO 110 I = 1, NNS
                  NRHS = NSVAL( I )
                  CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
                  CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                  IC = 0
                  S1 = SECOND( )
   90             CONTINUE
                  CALL CGTSV( M, NRHS, A, A( N ), A( 2*N ), B, LDB,
     $                        INFO )
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
                     CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 90
                  END IF
*
*                 Subtract the time used in CTIMMG.
*
                  ICL = 1
                  S1 = SECOND( )
  100             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
                     CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 100
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
                  OPS = SOPGB( 'CGTSV ', M, NRHS, 0, 0, IWORK )
                  RESLTS( I, IM, ILDA, 3 ) = SMFLOP( OPS, TIME, INFO )
  110          CONTINUE
  120       CONTINUE
         END IF
*
         IF( TIMSUB( 4 ) ) THEN
            CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
            CALL CTIMMG( 0, M, 1, B, N, 0, 0 )
            IC = 0
            S1 = SECOND( )
  130       CONTINUE
            CALL CGTSL( M, A, A( N ), A( 2*N ), B, INFO )
            S2 = SECOND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
               CALL CTIMMG( 0, M, 1, B, LDB, 0, 0 )
               GO TO 130
            END IF
*
*           Subtract the time used in CTIMMG.
*
            ICL = 1
            S1 = SECOND( )
  140       CONTINUE
            S2 = SECOND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               CALL CTIMMG( 14, M, M, A, 3*N, 0, 0 )
               CALL CTIMMG( 0, M, 1, B, LDB, 0, 0 )
               GO TO 140
            END IF
*
            TIME = ( TIME-UNTIME ) / REAL( IC )
            OPS = SOPGB( 'CGTSV ', M, 1, 0, 0, IWORK )
            RESLTS( 1, IM, 1, 4 ) = SMFLOP( OPS, TIME, INFO )
         END IF
  150 CONTINUE
*
*     Print a table of results for each timed routine.
*
      DO 170 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 170
         WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 .AND. ( TIMSUB( 2 ) .OR. TIMSUB( 3 ) ) ) THEN
            DO 160 I = 1, NLDA
               WRITE( NOUT, FMT = 9996 )I, LDAVAL( I )
  160       CONTINUE
         END IF
         WRITE( NOUT, FMT = * )
         IF( ISUB.EQ.1 ) THEN
            CALL SPRTBL( ' ', 'N', 1, LAVAL, NM, MVAL, 1, RESLTS, LDR1,
     $                   LDR2, NOUT )
         ELSE IF( ISUB.EQ.2 ) THEN
            WRITE( NOUT, FMT = 9998 )'N'
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 2 ), LDR1, LDR2, NOUT )
            WRITE( NOUT, FMT = 9998 )'T'
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 5 ), LDR1, LDR2, NOUT )
            WRITE( NOUT, FMT = 9998 )'C'
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 6 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.3 ) THEN
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 3 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.4 ) THEN
            CALL SPRTBL( ' ', 'N', 1, LAVAL, NM, MVAL, 1,
     $                   RESLTS( 1, 1, 1, 4 ), LDR1, LDR2, NOUT )
         END IF
  170 CONTINUE
*
  180 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( ' CGTTRS with TRANS = ''', A1, '''', / )
 9997 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9996 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of CTIMGT
*
      END
