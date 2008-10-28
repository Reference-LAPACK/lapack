      SUBROUTINE ZTIMTP( LINE, NN, NVAL, NNS, NSVAL, LA, TIMMIN, A, B,
     $                   RESLTS, LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LA, LDR1, LDR2, LDR3, NN, NNS, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            NSVAL( * ), NVAL( * )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, LDR3, * )
      COMPLEX*16         A( * ), B( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTIMTP times ZTPTRI and -TRS.
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
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix size N.
*
*  NNS     (input) INTEGER
*          The number of values of NRHS contained in the vector NSVAL.
*
*  NSVAL   (input) INTEGER array, dimension (NNS)
*          The values of the number of right hand sides NRHS.
*
*  LA      (input) INTEGER
*          The size of the arrays A and B.
*
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) COMPLEX*16 array, dimension (LA)
*
*  B       (workspace) COMPLEX*16 array, dimension (NMAX*NMAX)
*          where NMAX is the maximum value of N in NVAL.
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
*          The second dimension of RESLTS.  LDR2 >= max(1,NN).
*
*  LDR3    (input) INTEGER
*          The third dimension of RESLTS.  LDR3 >= 2.
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 2 )
*     ..
*     .. Local Scalars ..
      CHARACTER          UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, IN, INFO, ISUB, IUPLO, LDA, LDB,
     $                   MAT, N, NRHS
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          UPLOS( 2 )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            IDUMMY( 1 ), LAVAL( 1 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      LOGICAL            LSAME
      DOUBLE PRECISION   DMFLOP, DOPLA, DSECND
      EXTERNAL           LSAME, DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DPRTBL, ZTIMMG, ZTPTRI, ZTPTRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MOD
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'ZTPTRI', 'ZTPTRS' /
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'TP'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 100
*
*     Check that N*(N+1)/2 <= LA for the input values.
*
      CNAME = LINE( 1: 6 )
      LAVAL( 1 ) = LA
      CALL ATIMCK( 4, CNAME, NN, NVAL, 1, LAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 100
      END IF
*
*     Do first for UPLO = 'U', then for UPLO = 'L'
*
      DO 70 IUPLO = 1, 2
         UPLO = UPLOS( IUPLO )
         IF( LSAME( UPLO, 'U' ) ) THEN
            MAT = 12
         ELSE
            MAT = -12
         END IF
*
*        Do for each value of N:
*
         DO 60 IN = 1, NN
            N = NVAL( IN )
            LDA = N*( N+1 ) / 2
            LDB = N
            IF( MOD( N, 2 ).EQ.0 )
     $         LDB = LDB + 1
*
*           Time ZTPTRI
*
            IF( TIMSUB( 1 ) ) THEN
               CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
               IC = 0
               S1 = DSECND( )
   10          CONTINUE
               CALL ZTPTRI( UPLO, 'Non-unit', N, A, INFO )
               S2 = DSECND( )
               TIME = S2 - S1
               IC = IC + 1
               IF( TIME.LT.TIMMIN ) THEN
                  CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
                  GO TO 10
               END IF
*
*              Subtract the time used in ZTIMMG.
*
               ICL = 1
               S1 = DSECND( )
   20          CONTINUE
               S2 = DSECND( )
               UNTIME = S2 - S1
               ICL = ICL + 1
               IF( ICL.LE.IC ) THEN
                  CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
                  GO TO 20
               END IF
*
               TIME = ( TIME-UNTIME ) / DBLE( IC )
               OPS = DOPLA( 'ZTPTRI', N, N, 0, 0, 0 )
               RESLTS( 1, IN, IUPLO, 1 ) = DMFLOP( OPS, TIME, INFO )
            ELSE
*
*              Generate a triangular matrix A.
*
               CALL ZTIMMG( MAT, N, N, A, LDA, 0, 0 )
            END IF
*
*           Time ZTPTRS
*
            IF( TIMSUB( 2 ) ) THEN
               DO 50 I = 1, NNS
                  NRHS = NSVAL( I )
                  CALL ZTIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                  IC = 0
                  S1 = DSECND( )
   30             CONTINUE
                  CALL ZTPTRS( UPLO, 'No transpose', 'Non-unit', N,
     $                         NRHS, A, B, LDB, INFO )
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL ZTIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                     GO TO 30
                  END IF
*
*                 Subtract the time used in ZTIMMG.
*
                  ICL = 1
                  S1 = DSECND( )
   40             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL ZTIMMG( 0, N, NRHS, B, LDB, 0, 0 )
                     GO TO 40
                  END IF
*
                  TIME = ( TIME-UNTIME ) / DBLE( IC )
                  OPS = DOPLA( 'ZTPTRS', N, NRHS, 0, 0, 0 )
                  RESLTS( I, IN, IUPLO, 2 ) = DMFLOP( OPS, TIME, INFO )
   50          CONTINUE
            END IF
   60    CONTINUE
   70 CONTINUE
*
*     Print a table of results.
*
      DO 90 ISUB = 1, NSUBS
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 90
         WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         DO 80 IUPLO = 1, 2
            WRITE( NOUT, FMT = 9997 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )),
     $     UPLOS( IUPLO )
            IF( ISUB.EQ.1 ) THEN
               CALL DPRTBL( ' ', 'N', 1, IDUMMY, NN, NVAL, 1,
     $                      RESLTS( 1, 1, IUPLO, 1 ), LDR1, LDR2, NOUT )
            ELSE IF( ISUB.EQ.2 ) THEN
               CALL DPRTBL( 'NRHS', 'N', NNS, NSVAL, NN, NVAL, 1,
     $                      RESLTS( 1, 1, IUPLO, 2 ), LDR1, LDR2, NOUT )
            END IF
   80    CONTINUE
   90 CONTINUE
*
  100 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***', / )
 9997 FORMAT( 5X, A, ' with UPLO = ''', A1, '''', / )
      RETURN
*
*     End of ZTIMTP
*
      END
