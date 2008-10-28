      SUBROUTINE CTIMQ3( LINE, NM, MVAL, NVAL, NNB, NBVAL, NXVAL, NLDA,
     $                   LDAVAL, TIMMIN, A, COPYA, TAU, WORK, RWORK,
     $                   IWORK, RESLTS, LDR1, LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     Rewritten for timing qp3 code.
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, NLDA, NM, NNB, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NVAL( * ), NXVAL( * )
      REAL               RESLTS( LDR1, LDR2, * ), RWORK( * )
      COMPLEX            A( * ), COPYA( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  CTIMQ3 times the Rank-Revealing QR factorization of a
*  COMPLEX general matrix.
*
*  Two matrix types may be used for timing.  The number of types is
*  set in the parameter NMODE and the matrix types are set in the vector
*  MODES, using the following key:
*     2.  BREAK1    D(1:N-1)=1 and D(N)=1.0/COND in CLATMS
*     3.  GEOM      D(I)=COND**(-(I-1)/(N-1)) in CLATMS
*  These numbers are chosen to correspond with the matrix types in the
*  test code.
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
*          The number of values of M and N contained in the vectors
*          MVAL and NVAL.  The matrix sizes are used in pairs (M,N).
*
*  MVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix row dimension M.
*
*  NVAL    (input) INTEGER array, dimension (NM)
*          The values of the matrix column dimension N.
*
*  NK      (input) INTEGER
*          The number of values of K in the vector KVAL.
*
*  KVAL    (input) INTEGER array, dimension (NK)
*          The values of the matrix dimension K, used in SORMQR.
*
*  NNB     (input) INTEGER
*          The number of values of NB and NX contained in the
*          vectors NBVAL and NXVAL.  The blocking parameters are used
*          in pairs (NB,NX).
*
*  NBVAL   (input) INTEGER array, dimension (NNB)
*          The values of the blocksize NB.
*
*  NXVAL   (input) INTEGER array, dimension (NNB)
*          The values of the crossover point NX.
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
*  A       (workspace) COMPLEX array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  COPYA   (workspace) COMPLEX array, dimension (LDAMAX*NMAX)
*
*  WORK    (workspace) COMPLEX array, dimension (3*max(MMAX,NMAX))
*
*  RWORK   (workspace) REAL array, dimension (2*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  RESLTS  (workspace) REAL array, dimension
*                      (LDR1,LDR2,NLDA)
*          The timing results for each subroutine over the relevant
*          values of MODE, (M,N), and LDA.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NM).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
*
*
      INTEGER            NSUBS, NMODE
      PARAMETER          ( NSUBS = 1, NMODE = 2 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IM, IMODE, INB, INFO, LDA,
     $                   LW, M, MINMN, MODE, N, NB, NX
      REAL               COND, DMAX, OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), MODES( NMODE )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SLAMCH, SMFLOP, SOPLA
      EXTERNAL           SECOND, SLAMCH, SMFLOP, SOPLA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, CGEQP3, CLACPY, CLATMS, ICOPY,
     $                   SPRTB4, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CGEQP3' /
      DATA               MODES / 2, 3 /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'QP'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( .NOT.TIMSUB( 1 ) .OR. INFO.NE.0 )
     $   GO TO 90
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9996 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 90
      END IF
*
*     Set the condition number and scaling factor for the matrices
*     to be generated.
*
      DMAX = ONE
      COND = ONE / SLAMCH( 'Precision' )
*
*     Do for each type of matrix:
*
      DO 80 IMODE = 1, NMODE
         MODE = MODES( IMODE )
*
*
*        *****************
*        * Timing xGEQP3 *
*        *****************
*
*        Do for each value of LDA:
*
         DO 60 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
*
*           Do for each pair of values (M,N):
*
            DO 50 IM = 1, NM
               M = MVAL( IM )
               N = NVAL( IM )
               MINMN = MIN( M, N )
*
*              Generate a test matrix of size m by n using the
*              singular value distribution indicated by MODE.
*
               CALL CLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', RWORK,
     $                      MODE, COND, DMAX, M, N, 'No packing', COPYA,
     $                      LDA, WORK, INFO )
*
*              Do for each pair of values ( NB, NX ) in NBVAL and NXVAL:
*
               DO 40 INB = 1, NNB
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
                  NX = NXVAL( INB )
                  CALL XLAENV( 3, NX )
*
*
*                 CGEQP3
*
                  LW = MAX( 1, ( N+1 )*NB )
                  DO 10 I = 1, N
                     IWORK( N+I ) = 0
   10             CONTINUE
                  IC = 0
                  S1 = SECOND( )
   20             CONTINUE
                  CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                  CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL CGEQP3( M, N, A, LDA, IWORK, TAU, WORK, LW,
     $                         RWORK, INFO )
*
                  IF( INFO.NE.0 ) THEN
                     WRITE( *, FMT = * )'>>>Warning: INFO returned by ',
     $                  'CGEQP3 is:', INFO
                     INFO = 0
                  END IF
*
                  S2 = SECOND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 20
*
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN ) THEN
                     CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                     CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                     GO TO 20
                  END IF
*
*                 Subtract the time used in CLACPY.
*
                  ICL = 1
                  S1 = SECOND( )
   30             CONTINUE
                  S2 = SECOND( )
                  UNTIME = S2 - S1
                  ICL = ICL + 1
                  IF( ICL.LE.IC ) THEN
                     CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                     CALL CLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                     GO TO 30
                  END IF
*
                  TIME = ( TIME-UNTIME ) / REAL( IC )
*
*                 The number of flops of yGEQP3 is approximately the
*                 the number of flops of yGEQPF.
*
                  OPS = SOPLA( 'CGEQPF', M, N, 0, 0, NB )
                  RESLTS( INB, IM, ILDA ) = SMFLOP( OPS, TIME, INFO )
*
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
*
*        Print the results for each value of K and type of matrix.
*
         WRITE( NOUT, FMT = 9999 )
     $     SUBNAM( 1 )(1:ILA_LEN_TRIM( SUBNAM( 1 ) ))
         WRITE( NOUT, FMT = 9998 )IMODE
         DO 70 I = 1, NLDA
            WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
   70    CONTINUE
         WRITE( NOUT, FMT = * )
         CALL SPRTB4( '(  NB,  NX)', 'M', 'N', NNB, NBVAL, NXVAL, NM,
     $                MVAL, NVAL, NLDA, RESLTS( 1, 1, 1 ), LDR1, LDR2,
     $                NOUT )
*
   80 CONTINUE
*
 9999 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9998 FORMAT( 5X, 'type of matrix:', I4 )
 9997 FORMAT( 5X, 'line ', I4, ' with LDA = ', I4 )
 9996 FORMAT( 1X, A, ' timing run not attempted', / )
*
   90 CONTINUE
      RETURN
*
*     End of CTIMQ3
*
      END
