      SUBROUTINE CTIMPT( LINE, NM, MVAL, NNS, NSVAL, NLDA, LDAVAL,
     $                   TIMMIN, D, E, B, RESLTS, LDR1, LDR2, LDR3,
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
      INTEGER            LDAVAL( * ), MVAL( * ), NSVAL( * )
      REAL               D( * ), RESLTS( LDR1, LDR2, LDR3, * )
      COMPLEX            B( * ), E( * )
*     ..
*
*  Purpose
*  =======
*
*  CTIMPT times CPTTRF, -TRS, -SV, and -SL.
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
*  D       (workspace) REAL array, dimension (NMAX)
*          where NMAX is the maximum value permitted for N.
*
*  E       (workspace) COMPLEX array, dimension (2*NMAX)
*          where NMAX is the maximum value permitted for N.
*
*  B       (workspace) COMPLEX array, dimension (LDAMAX*NMAX)
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
      CHARACTER          UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IM, INFO, IOFF, ISUB, IUPLO,
     $                   J, LDB, M, N, NRHS
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), LAVAL( 1 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SMFLOP, SOPLA
      EXTERNAL           SECOND, SMFLOP, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, CLARNV, CPTSL, CPTSV, CPTTRF,
     $                   CPTTRS, CTIMMG, SPRTBL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CPTTRF', 'CPTTRS', 'CPTSV ',
     $                   'CPTSL ' /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'PT'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 280
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
      DO 250 IM = 1, NM
*
         M = MVAL( IM )
         N = MAX( M, 1 )
*
*        Time CPTTRF
*
         IF( TIMSUB( 1 ) ) THEN
            DO 20 J = 1, N
               D( J ) = 3.0
   20       CONTINUE
            CALL CLARNV( 2, ISEED, N-1, E )
            IC = 0
            S1 = SECOND( )
   30       CONTINUE
            CALL CPTTRF( M, D, E, INFO )
            S2 = SECOND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               DO 40 J = 1, N
                  D( J ) = 3.0
   40          CONTINUE
               CALL CLARNV( 2, ISEED, N-1, E )
               GO TO 30
            END IF
*
*           Subtract the time used in STIMMG.
*
            ICL = 1
            S1 = SECOND( )
   50       CONTINUE
            S2 = SECOND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               DO 60 J = 1, N
                  D( J ) = 3.0
   60          CONTINUE
               CALL CLARNV( 2, ISEED, N-1, E )
               GO TO 50
            END IF
*
            TIME = ( TIME-UNTIME ) / REAL( IC )
            OPS = SOPLA( 'CPTTRF', M, 0, 0, 0, 0 )
            RESLTS( 1, IM, 1, 1 ) = SMFLOP( OPS, TIME, INFO )
*
         ELSE
            IC = 0
            DO 70 J = 1, N
               D( J ) = 3.0
   70       CONTINUE
            CALL CLARNV( 2, ISEED, N-1, E )
         END IF
*
*        Generate another matrix and factor it using CPTTRF so
*        that the factored form can be used in timing the other
*        routines.
*
         IF( IC.NE.1 )
     $      CALL CPTTRF( M, D, E, INFO )
*
*        Time CPTTRS
*
         DO 190 IUPLO = 1, 2
            IF( IUPLO.EQ.1 ) THEN
               UPLO = 'U'
               IOFF = 0
            ELSE
               UPLO = 'L'
               IOFF = 3
            END IF
            IF( TIMSUB( 2 ) ) THEN
               DO 110 ILDA = 1, NLDA
                  LDB = LDAVAL( ILDA )
                  DO 100 I = 1, NNS
                     NRHS = NSVAL( I )
                     CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
   80                CONTINUE
                     CALL CPTTRS( UPLO, M, NRHS, D, E, B, LDB, INFO )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                        GO TO 80
                     END IF
*
*                    Subtract the time used in STIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
   90                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                        GO TO 90
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPLA( 'CPTTRS', M, NRHS, 0, 0, 0 )
                     RESLTS( I, IM, ILDA, IOFF+2 ) = SMFLOP( OPS, TIME,
     $                  INFO )
  100             CONTINUE
  110          CONTINUE
            END IF
*
            IF( TIMSUB( 3 ) .AND. IUPLO.EQ.1 ) THEN
               DO 180 ILDA = 1, NLDA
                  LDB = LDAVAL( ILDA )
                  DO 170 I = 1, NNS
                     NRHS = NSVAL( I )
                     DO 120 J = 1, N
                        D( J ) = 3.0
  120                CONTINUE
                     CALL CLARNV( 2, ISEED, N-1, E )
                     CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     IC = 0
                     S1 = SECOND( )
  130                CONTINUE
*                     CALL CPTSV( UPLO, M, NRHS, D, E, B, LDB, INFO )
                     CALL CPTSV( M, NRHS, D, E, B, LDB, INFO )
                     S2 = SECOND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        DO 140 J = 1, N
                           D( J ) = 3.0
  140                   CONTINUE
                        CALL CLARNV( 2, ISEED, N-1, E )
                        CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                        GO TO 130
                     END IF
*
*                    Subtract the time used in CTIMMG.
*
                     ICL = 1
                     S1 = SECOND( )
  150                CONTINUE
                     S2 = SECOND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        DO 160 J = 1, N
                           D( J ) = 3.0
  160                   CONTINUE
                        CALL CLARNV( 2, ISEED, N-1, E )
                        CALL CTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                        GO TO 150
                     END IF
*
                     TIME = ( TIME-UNTIME ) / REAL( IC )
                     OPS = SOPLA( 'CPTSV ', M, NRHS, 0, 0, 0 )
                     RESLTS( I, IM, ILDA, IOFF+3 ) = SMFLOP( OPS, TIME,
     $                  INFO )
  170             CONTINUE
  180          CONTINUE
            END IF
  190    CONTINUE
*
         IF( TIMSUB( 4 ) ) THEN
            DO 200 J = 1, N
               E( J ) = 3.0
  200       CONTINUE
            CALL CLARNV( 2, ISEED, N-1, E( N+1 ) )
            CALL CTIMMG( 0, M, 1, B, N, 0, 0 )
            IC = 0
            S1 = SECOND( )
  210       CONTINUE
            CALL CPTSL( M, E, E( N+1 ), B )
            S2 = SECOND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               DO 220 J = 1, N
                  E( J ) = 3.0
  220          CONTINUE
               CALL CLARNV( 2, ISEED, N-1, E( N+1 ) )
               CALL CTIMMG( 0, M, 1, B, N, 0, 0 )
               GO TO 210
            END IF
*
*           Subtract the time used in STIMMG.
*
            ICL = 1
            S1 = SECOND( )
  230       CONTINUE
            S2 = SECOND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               DO 240 J = 1, N
                  E( J ) = 3.0
  240          CONTINUE
               CALL CLARNV( 2, ISEED, N-1, E( N+1 ) )
               CALL CTIMMG( 0, M, 1, B, N, 0, 0 )
               GO TO 230
            END IF
*
            TIME = ( TIME-UNTIME ) / REAL( IC )
            OPS = SOPLA( 'CPTSV ', M, 1, 0, 0, 0 )
            RESLTS( 1, IM, 1, 4 ) = SMFLOP( OPS, TIME, INFO )
         END IF
  250 CONTINUE
*
*     Print a table of results for each timed routine.
*
      DO 270 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 270
         WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 .AND. ( TIMSUB( 2 ) .OR. TIMSUB( 3 ) ) ) THEN
            DO 260 I = 1, NLDA
               WRITE( NOUT, FMT = 9996 )I, LDAVAL( I )
  260       CONTINUE
         END IF
         WRITE( NOUT, FMT = * )
         IF( ISUB.EQ.1 ) THEN
            CALL SPRTBL( ' ', 'N', 1, LAVAL, NM, MVAL, 1, RESLTS, LDR1,
     $                   LDR2, NOUT )
         ELSE IF( ISUB.EQ.2 ) THEN
            WRITE( NOUT, FMT = 9998 )'CPTTRS', 'U'
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 2 ), LDR1, LDR2, NOUT )
            WRITE( NOUT, FMT = 9998 )'CPTTRS', 'L'
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 5 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.3 ) THEN
*            WRITE( NOUT, FMT = 126 ) 'CPTSV ', 'U'
            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 3 ), LDR1, LDR2, NOUT )
*            WRITE( NOUT, FMT = 126 ) 'CPTSV ', 'L'
*            CALL SPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
*     $                   RESLTS( 1, 1, 1, 6 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.4 ) THEN
            CALL SPRTBL( ' ', 'N', 1, LAVAL, NM, MVAL, 1,
     $                   RESLTS( 1, 1, 1, 4 ), LDR1, LDR2, NOUT )
         END IF
  270 CONTINUE
*
  280 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( 1X, A, ' with UPLO = ''', A1, '''', / )
 9997 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9996 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of CTIMPT
*
      END
