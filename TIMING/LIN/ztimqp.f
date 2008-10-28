      SUBROUTINE ZTIMQP( LINE, NM, MVAL, NVAL, NLDA, LDAVAL, TIMMIN, A,
     $                   COPYA, TAU, WORK, RWORK, IWORK, RESLTS, LDR1,
     $                   LDR2, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, NLDA, NM, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NVAL( * )
      DOUBLE PRECISION   RESLTS( LDR1, LDR2, * ), RWORK( * )
      COMPLEX*16         A( * ), COPYA( * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZTIMQP times the LAPACK routines to perform the QR factorization with
*  column pivoting of a COMPLEX*16 general matrix.
*
*  Two matrix types may be used for timing.  The number of types is
*  set in the parameter NMODE and the matrix types are set in the vector
*  MODES, using the following key:
*     2.  BREAK1    D(1:N-1)=1 and D(N)=1.0/COND in ZLATMS
*     3.  GEOM      D(I)=COND**(-(I-1)/(N-1)) in ZLATMS
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
*  NLDA    (input) INTEGER
*          The number of values of LDA contained in the vector LDAVAL.
*
*  LDAVAL  (input) INTEGER array, dimension (NLDA)
*          The values of the leading dimension of the array A.
*
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  COPYA   (workspace) COMPLEX*16 array, dimension (LDAMAX*NMAX)
*
*  TAU     (workspace) COMPLEX*16 array, dimension (min(M,N))
*
*  WORK    (workspace) COMPLEX*16 array, dimension (3*max(MMAX,NMAX))
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*NMAX)
*
*  IWORK   (workspace) INTEGER array, dimension (2*NMAX)
*
*  RESLTS  (workspace) DOUBLE PRECISION array, dimension
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
      INTEGER            NSUBS, NMODE
      PARAMETER          ( NSUBS = 1, NMODE = 2 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, IC, ICL, ILDA, IM, IMODE, INFO, LDA, M,
     $                   MINMN, MODE, N
      DOUBLE PRECISION   COND, DMAX, OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), MODES( NMODE )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      DOUBLE PRECISION   DLAMCH, DMFLOP, DOPLA, DSECND
      EXTERNAL           DLAMCH, DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DPRTB5, ICOPY, ZGEQPF, ZLACPY,
     $                   ZLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'ZGEQPF' /
      DATA               MODES / 2, 3 /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'QP'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( .NOT.TIMSUB( 1 ) .OR. INFO.NE.0 )
     $   GO TO 80
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 1, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 80
      END IF
*
*     Set the condition number and scaling factor for the matrices
*     to be generated.
*
      DMAX = ONE
      COND = ONE / DLAMCH( 'Precision' )
*
*     Do for each pair of values (M,N):
*
      DO 60 IM = 1, NM
         M = MVAL( IM )
         N = NVAL( IM )
         MINMN = MIN( M, N )
*
*        Do for each value of LDA:
*
         DO 50 ILDA = 1, NLDA
            LDA = LDAVAL( ILDA )
            DO 40 IMODE = 1, NMODE
               MODE = MODES( IMODE )
*
*              Generate a test matrix of size m by n using the
*              singular value distribution indicated by MODE.
*
               DO 10 I = 1, N
                  IWORK( N+I ) = 0
   10          CONTINUE
               CALL ZLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', RWORK,
     $                      MODE, COND, DMAX, M, N, 'No packing', COPYA,
     $                      LDA, WORK, INFO )
*
*              ZGEQPF:  QR factorization with column pivoting
*
               CALL ZLACPY( 'All', M, N, COPYA, LDA, A, LDA )
               CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
               IC = 0
               S1 = DSECND( )
   20          CONTINUE
               CALL ZGEQPF( M, N, A, LDA, IWORK, TAU, WORK, RWORK,
     $                      INFO )
               S2 = DSECND( )
               TIME = S2 - S1
               IC = IC + 1
               IF( TIME.LT.TIMMIN ) THEN
                  CALL ZLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                  GO TO 20
               END IF
*
*              Subtract the time used in ZLACPY and ICOPY.
*
               ICL = 1
               S1 = DSECND( )
   30          CONTINUE
               S2 = DSECND( )
               UNTIME = S2 - S1
               ICL = ICL + 1
               IF( ICL.LE.IC ) THEN
                  CALL ZLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL ICOPY( N, IWORK( N+1 ), 1, IWORK, 1 )
                  GO TO 30
               END IF
*
               TIME = ( TIME-UNTIME ) / DBLE( IC )
               OPS = DOPLA( 'ZGEQPF', M, N, 0, 0, 1 )
               RESLTS( IMODE, IM, ILDA ) = DMFLOP( OPS, TIME, INFO )
*
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
*
*     Print tables of results
*
      WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( 1 )(1:ILA_LEN_TRIM( SUBNAM( 1 ) ))
      IF( NLDA.GT.1 ) THEN
         DO 70 I = 1, NLDA
            WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
   70    CONTINUE
      END IF
      WRITE( NOUT, FMT = * )
      CALL DPRTB5( 'Type', 'M', 'N', NMODE, MODES, NM, MVAL, NVAL, NLDA,
     $             RESLTS, LDR1, LDR2, NOUT )
   80 CONTINUE
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops ***' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
      RETURN
*
*     End of ZTIMQP
*
      END
