      SUBROUTINE STIMMM( VNAME, LAB2, NN, NVAL, NLDA, LDAVAL, TIMMIN, A,
     $                   B, C, RESLTS, LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    LAB2, VNAME
      INTEGER            LDR1, LDR2, NLDA, NN, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            LDAVAL( * ), NVAL( * )
      REAL               A( * ), B( * ), C( * ), RESLTS( LDR1, LDR2, * )
*     ..
*
*  Purpose
*  =======
*
*  STIMMM times SGEMM.
*
*  Arguments
*  =========
*
*  VNAME   (input) CHARACTER*(*)
*          The name of the Level 3 BLAS routine to be timed.
*
*  LAB2    (input) CHARACTER*(*)
*          The name of the variable given in NVAL.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix dimension N.
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
*  A       (workspace) REAL array, dimension (LDAMAX*NMAX)
*             where LDAMAX and NMAX are the maximum values permitted
*             for LDA and N.
*
*  B       (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  C       (workspace) REAL array, dimension (LDAMAX*NMAX)
*
*  RESLTS  (output) REAL array, dimension (LDR1,LDR2,NLDA)
*          The timing results for each subroutine over the relevant
*          values of N and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= 1.
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NN).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      REAL               ONE
      PARAMETER          ( NSUBS = 1, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IN, INFO, ISUB, LDA, N
      REAL               OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            IDUMMY( 1 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      LOGICAL            LSAMEN
      REAL               SECOND, SMFLOP, SOPBL3
      EXTERNAL           LSAMEN, SECOND, SMFLOP, SOPBL3
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, SGEMM, SPRTBL, STIMMG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'SGEMM ' /
*     ..
*     .. Executable Statements ..
*
      CNAME = VNAME
      DO 10 ISUB = 1, NSUBS
         TIMSUB( ISUB ) = LSAMEN( 6, CNAME, SUBNAM( ISUB ) )
         IF( TIMSUB( ISUB ) )
     $      GO TO 20
   10 CONTINUE
      WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
      GO TO 80
   20 CONTINUE
*
*     Check that N <= LDA for the input values.
*
      CALL ATIMCK( 2, CNAME, NN, NVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9998 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 80
      END IF
*
      DO 60 ILDA = 1, NLDA
         LDA = LDAVAL( ILDA )
         DO 50 IN = 1, NN
            N = NVAL( IN )
*
*           Time SGEMM
*
            CALL STIMMG( 1, N, N, A, LDA, 0, 0 )
            CALL STIMMG( 0, N, N, B, LDA, 0, 0 )
            CALL STIMMG( 1, N, N, C, LDA, 0, 0 )
            IC = 0
            S1 = SECOND( )
   30       CONTINUE
            CALL SGEMM( 'No transpose', 'No transpose', N, N, N, ONE, A,
     $                  LDA, B, LDA, ONE, C, LDA )
            S2 = SECOND( )
            TIME = S2 - S1
            IC = IC + 1
            IF( TIME.LT.TIMMIN ) THEN
               CALL STIMMG( 1, N, N, C, LDA, 0, 0 )
               GO TO 30
            END IF
*
*           Subtract the time used in STIMMG.
*
            ICL = 1
            S1 = SECOND( )
   40       CONTINUE
            S2 = SECOND( )
            UNTIME = S2 - S1
            ICL = ICL + 1
            IF( ICL.LE.IC ) THEN
               CALL STIMMG( 1, N, N, C, LDA, 0, 0 )
               GO TO 40
            END IF
*
            TIME = ( TIME-UNTIME ) / REAL( IC )
            OPS = SOPBL3( 'SGEMM ', N, N, N )
            RESLTS( 1, IN, ILDA ) = SMFLOP( OPS, TIME, 0 )
   50    CONTINUE
   60 CONTINUE
*
*     Print the table of results on unit NOUT.
*
      WRITE( NOUT, FMT = 9997 )VNAME
      IF( NLDA.EQ.1 ) THEN
         WRITE( NOUT, FMT = 9996 )LDAVAL( 1 )
      ELSE
         DO 70 I = 1, NLDA
            WRITE( NOUT, FMT = 9995 )I, LDAVAL( I )
   70    CONTINUE
      END IF
      WRITE( NOUT, FMT = * )
      CALL SPRTBL( ' ', LAB2, 1, IDUMMY, NN, NVAL, NLDA, RESLTS, LDR1,
     $             LDR2, NOUT )
*
   80 CONTINUE
      RETURN
 9999 FORMAT( 1X, A, ':  Unrecognized path or subroutine name', / )
 9998 FORMAT( 1X, A, ' timing run not attempted', / )
 9997 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9996 FORMAT( 5X, 'with LDA = ', I5 )
 9995 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
*
*     End of STIMMM
*
      END
