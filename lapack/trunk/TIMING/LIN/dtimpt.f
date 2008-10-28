      SUBROUTINE DTIMPT( LINE, NM, MVAL, NNS, NSVAL, NLDA, LDAVAL,
     $                   TIMMIN, A, B, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NLDA, NM, NNS, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            LDAVAL( * ), MVAL( * ), NSVAL( * )
      DOUBLE PRECISION   A( * ), B( * ), RESLTS( LDR1, LDR2, LDR3, * )
*     ..
*
*  Purpose
*  =======
*
*  DTIMPT times DPTTRF, -TRS, -SV, and -SL.
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
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) DOUBLE PRECISION array, dimension (NMAX*2)
*          where NMAX is the maximum value permitted for N.
*
*  B       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*
*  RESLTS  (output) DOUBLE PRECISION array, dimension
*                   (LDR1,LDR2,LDR3,NSUBS)
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
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IM, INFO, ISUB, LDB, M, N,
     $                   NRHS
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            LAVAL( 1 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      DOUBLE PRECISION   DMFLOP, DOPLA, DSECND
      EXTERNAL           DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DPRTBL, DPTSV, DPTTRF, DPTTRS,
     $                   DTIMMG, DPTSL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'DPTTRF', 'DPTTRS', 'DPTSV ',
     $                   'DPTSL ' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'PT'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 170
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
      DO 140 IM = 1, NM
*
         M = MVAL( IM )
         N = MAX( M, 1 )
*
*        Time DPTTRF
*
         IF( TIMSUB( 1 ) ) THEN
            CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
            IC = 0
            S1 = DSECND( )
   20       CONTINUE
            CALL DPTTRF( M, A, A( N+1 ), INFO )
            S2 = DSECND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
               GO TO 20
            END IF
*
*           Subtract the time used in DTIMMG.
*
            ICL = 1
            S1 = DSECND( )
   30       CONTINUE
            S2 = DSECND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
               GO TO 30
            END IF
*
            TIME = ( TIME-UNTIME ) / DBLE( IC )
            OPS = DOPLA( 'DPTTRF', M, 0, 0, 0, 0 )
            RESLTS( 1, IM, 1, 1 ) = DMFLOP( OPS, TIME, INFO )
*
         ELSE
            IC = 0
            CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
         END IF
*
*        Generate another matrix and factor it using DPTTRF so
*        that the factored form can be used in timing the other
*        routines.
*
         IF( IC.NE.1 )
     $      CALL DPTTRF( M, A, A( N+1 ), INFO )
*
*        Time DPTTRS
*
         IF( TIMSUB( 2 ) ) THEN
            DO 70 ILDA = 1, NLDA
               LDB = LDAVAL( ILDA )
               DO 60 I = 1, NNS
                  NRHS = NSVAL( I )
                  CALL DTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                  IC = 0
                  S1 = DSECND( )
   40             CONTINUE
                  CALL DPTTRS( M, NRHS, A, A( N+1 ), B, LDB, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL DTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 40
                  END IF
*
*                 Subtract the time used in DTIMMG.
*
                  ICL = 1
                  S1 = DSECND( )
   50             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL DTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 50
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'DPTTRS', M, NRHS, 0, 0, 0 )
                  RESLTS( I, IM, ILDA, 2 ) = DMFLOP( OPS, TIME, INFO )
   60          CONTINUE
   70       CONTINUE
         END IF
*
         IF( TIMSUB( 3 ) ) THEN
            DO 110 ILDA = 1, NLDA
               LDB = LDAVAL( ILDA )
               DO 100 I = 1, NNS
                  NRHS = NSVAL( I )
                  CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
                  CALL DTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                  IC = 0
                  S1 = DSECND( )
   80             CONTINUE
                  CALL DPTSV( M, NRHS, A, A( N+1 ), B, LDB, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
                     CALL DTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 80
                  END IF
*
*                 Subtract the time used in DTIMMG.
*
                  ICL = 1
                  S1 = DSECND( )
   90             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
                     CALL DTIMMG( 0, M, NRHS, B, LDB, 0, 0 )
                     GO TO 90
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'DPTSV ', M, NRHS, 0, 0, 0 )
                  RESLTS( I, IM, ILDA, 3 ) = DMFLOP( OPS, TIME, INFO )
  100          CONTINUE
  110       CONTINUE
         END IF
*
         IF( TIMSUB( 4 ) ) THEN
            CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
            CALL DTIMMG( 0, M, 1, B, N, 0, 0 )
            IC = 0
            S1 = DSECND( )
  120       CONTINUE
            CALL DPTSL( M, A, A( N+1 ), B )
            S2 = DSECND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
               CALL DTIMMG( 0, M, 1, B, N, 0, 0 )
               GO TO 120
            END IF
*
*           Subtract the time used in DTIMMG.
*
            ICL = 1
            S1 = DSECND( )
  130       CONTINUE
            S2 = DSECND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               CALL DTIMMG( 13, M, M, A, 2*N, 0, 0 )
               CALL DTIMMG( 0, M, 1, B, N, 0, 0 )
               GO TO 130
            END IF
*
            TIME = ( TIME-UNTIME ) / DBLE( IC )
            OPS = DOPLA( 'DPTSV ', M, 1, 0, 0, 0 )
            RESLTS( 1, IM, 1, 4 ) = DMFLOP( OPS, TIME, INFO )
         END IF
  140 CONTINUE
*
*     Print a table of results for each timed routine.
*
      DO 160 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 160
         WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 .AND. ( TIMSUB( 2 ) .OR. TIMSUB( 3 ) ) ) THEN
            DO 150 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  150       CONTINUE
         END IF
         WRITE( NOUT, FMT = * )
         IF( ISUB.EQ.1 ) THEN
            CALL DPRTBL( ' ', 'N', 1, LAVAL, NM, MVAL, 1, RESLTS, LDR1,
     $                   LDR2, NOUT )
         ELSE IF( ISUB.EQ.2 ) THEN
            CALL DPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 2 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.3 ) THEN
            CALL DPRTBL( 'NRHS', 'N', NNS, NSVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, 3 ), LDR1, LDR2, NOUT )
         ELSE IF( ISUB.EQ.4 ) THEN
            CALL DPRTBL( ' ', 'N', 1, LAVAL, NM, MVAL, 1,
     $                   RESLTS( 1, 1, 1, 4 ), LDR1, LDR2, NOUT )
         END IF
  160 CONTINUE
*
  170 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of DTIMPT
*
      END
