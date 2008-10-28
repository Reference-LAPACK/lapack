      SUBROUTINE CTIMLS( LINE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB,
     $                   NBVAL, NXVAL, NLDA, LDAVAL, TIMMIN, A, COPYA,
     $                   B, COPYB, S, COPYS, OPCTBL, TIMTBL, FLPTBL,
     $                   WORK, RWORK, IWORK, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            NLDA, NM, NN, NNB, NNS, NOUT
      REAL               TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            IWORK( * ), LDAVAL( * ), MVAL( * ), NBVAL( * ),
     $                   NSVAL( * ), NVAL( * ), NXVAL( * )
      REAL               COPYS( * ),
     $                   FLPTBL( 6, 6, NM*NN*NNS*NLDA*(NNB+1), * ),
     $                   OPCTBL( 6, 6, NM*NN*NNS*NLDA*(NNB+1), * ),
     $                   S( * ), RWORK( * ),
     $                   TIMTBL( 6, 6, NM*NN*NNS*NLDA*(NNB+1), * ) 
      COMPLEX            A( * ), B( * ), COPYA( * ), COPYB( * ),
     $                   WORK( * )
*     ..
*     .. Common blocks ..
      COMMON             / LSTIME / OPCNT, TIMNG
*     ..
*     .. Arrays in Common ..
      REAL               OPCNT( 6 ), TIMNG( 6 )
*     ..
*
*  Purpose
*  =======
*
*  CTIMLS times the least squares driver routines CGELS, CGELSS, CGELSX,
*  CGELSY and CGELSD.
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
*          The values of the matrix row dimension M.
*
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
*
*  NNS     (input) INTEGER
*          The number of values of NRHS contained in the vector NSVAL.
*
*  NSVAL   (input) INTEGER array, dimension (NNS)
*          The values of the number of right hand sides NRHS.
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
*  A       (workspace) COMPLEX array, dimension (MMAX*NMAX)
*          where MMAX is the maximum value of M in MVAL and NMAX is the
*          maximum value of N in NVAL.
*
*  COPYA   (workspace) COMPLEX array, dimension (MMAX*NMAX)
*
*  B       (workspace) COMPLEX array, dimension (MMAX*NSMAX)
*          where MMAX is the maximum value of M in MVAL and NSMAX is the
*          maximum value of NRHS in NSVAL.
*
*  COPYB   (workspace) COMPLEX array, dimension (MMAX*NSMAX)
*
*  S       (workspace) REAL array, dimension
*                      (min(MMAX,NMAX))
*
*  COPYS   (workspace) REAL array, dimension
*                      (min(MMAX,NMAX))
*
*  OPCTBL  (workspace) REAL array, dimension
*                      (6,6,(NNB+1)*NLDA,NM*NN*NNS,5)
*
*  TIMTBL  (workspace) REAL array, dimension
*                      (6,6,(NNB+1)*NLDA,NM*NN*NNS,5)
*
*  FLPTBL  (workspace) REAL array, dimension
*                      (6,6,(NNB+1)*NLDA,NM*NN*NNS,5)
*
*  WORK    (workspace) COMPLEX array,
*                      dimension (MMAX*NMAX + 4*NMAX + MMAX).
*
*  IWORK   (workspace) INTEGER array, dimension (NMAX)
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MTYPE, NSUBS
      PARAMETER          ( MTYPE = 6, NSUBS = 5 )
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E0, ZERO = 0.0E0 )
      COMPLEX            CONE, CZERO
      PARAMETER          ( CONE = ( 1.0E0, 0.0E0 ),
     $                     CZERO = ( 0.0E0, 0.0E0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANS
      CHARACTER*3        PATH
      INTEGER            CRANK, I, ILDA, IM, IN, INB, INFO, INS, IRANK,
     $                   ISCALE, ISUB, ITBL, ITRAN, ITYPE, LDA, LDB,
     $                   LDWORK, LWLSY, LWORK, M, MNMIN, N, NB,
     $                   NCLS, NCLSD, NCLSS, NCLSX, NCLSY,
     $                   NCALL, NCOLS, NRHS, NROWS, RANK
      REAL               EPS, NORMA, NORMB, RCOND, S1, S2, TIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), ISEEDY( 4 ), NDATA( NSUBS )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      REAL               SECOND, SASUM, SLAMCH, SMFLOP
      EXTERNAL           SECOND, SASUM, SLAMCH, SMFLOP
*     ..
*     .. External Subroutines ..
      EXTERNAL           CSSCAL, CGELS, CGELSD,
     $                   CGELSS, CGELSX, CGELSY, CGEMM,
     $                   CLACPY, CLARNV, CQRT13, CQRT15,
     $                   SCOPY, SLASET, SPRTLS, SSCAL,
     $                   XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          REAL, MAX, MIN, SQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER(32)      SRNAMT
      INTEGER            INFOT, IOUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'CGELS ', 'CGELSX', 'CGELSY',
     $                            'CGELSS', 'CGELSD' /
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               NDATA  / 4, 6, 6, 6, 5 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'LS'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 230
*
*     Initialize constants and the random number seed.
*
      NCLS = 0
      NCLSD = 0
      NCLSS = 0
      NCLSX = 0
      NCLSY = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = SLAMCH( 'Epsilon' )
*
*     Threshold for rank estimation
*
      RCOND = SQRT( EPS ) - ( SQRT( EPS )-EPS ) / 2
*
      INFOT = 0
      CALL XLAENV( 2, 2 )
*
      DO 200 IM = 1, NM
         M = MVAL( IM )
*
         DO 190 IN = 1, NN
            N = NVAL( IN )
            MNMIN = MIN( M, N )
*
            DO 180 INS = 1, NNS
               NRHS = NSVAL( INS )
               LWORK = MAX( 1, ( M+NRHS )*( N+2 ), ( N+NRHS )*( M+2 ),
     $                 M*N+4*MNMIN+MAX( M, N ), 2*N+M )
*
               DO 170 ILDA = 1, NLDA
                  LDA = MAX( 1, LDAVAL( ILDA ) )
                  LDB = MAX( 1, LDAVAL( ILDA ), M, N )
*
                  DO 160 IRANK = 1, 2
*
                     DO 150 ISCALE = 1, 3
*
                        IF( IRANK.EQ.1 .AND. TIMSUB( 1 ) ) THEN
*
*                          Time CGELS
*
*                          Generate a matrix of scaling type ISCALE
*
                           CALL CQRT13( ISCALE, M, N, COPYA, LDA,
     $                               NORMA, ISEED )
                           DO 50 INB = 1, NNB
                              NB = NBVAL( INB )
                              CALL XLAENV( 1, NB )
                              CALL XLAENV( 3, NXVAL( INB ) )
*
                              DO 40 ITRAN = 1, 2
                                 ITYPE = ( ITRAN-1 )*3 + ISCALE 
                                 IF( ITRAN.EQ.1 ) THEN
                                    TRANS = 'N'
                                    NROWS = M
                                    NCOLS = N
                                 ELSE
                                    TRANS = 'C'
                                    NROWS = N
                                    NCOLS = M
                                 END IF
                                 LDWORK = MAX( 1, NCOLS )
*
*                                Set up a consistent rhs
*
                                 IF( NCOLS.GT.0 ) THEN
                                    CALL CLARNV( 2, ISEED, NCOLS*NRHS,
     $                                           WORK )
                                    CALL CSSCAL( NCOLS*NRHS,
     $                                           ONE / REAL( NCOLS ),
     $                                           WORK, 1 )
                                 END IF
                                 CALL CGEMM( TRANS, 'No transpose',
     $                                       NROWS, NRHS, NCOLS, CONE,
     $                                       COPYA, LDA, WORK, LDWORK,
     $                                       CZERO, B, LDB )
                                 CALL CLACPY( 'Full', NROWS, NRHS, B,
     $                                        LDB, COPYB, LDB )
*
*                                Solve LS or overdetermined system
*
                                 NCALL = 0
                                 TIME = ZERO
                                 CALL SLASET( 'Full', NDATA( 1 ), 1,
     $                                        ZERO, ZERO, OPCNT,
     $                                        NDATA( 1 ) )
                                 CALL SLASET( 'Full', NDATA( 1 ), 1,
     $                                        ZERO, ZERO, TIMNG,
     $                                        NDATA( 1 ) )
   20                            CONTINUE
                                 IF( M.GT.0 .AND. N.GT.0 ) THEN
                                 CALL CLACPY( 'Full', M, N, COPYA, LDA,
     $                                        A, LDA )
                                 CALL CLACPY( 'Full', NROWS, NRHS,
     $                                        COPYB, LDB, B, LDB )
                                 END IF
                                 SRNAMT = 'CGELS '
                                 NCALL = NCALL + 1
                                 S1 = SECOND( )
                                 CALL CGELS( TRANS, M, N, NRHS, A, LDA,
     $                                       B, LDB, WORK, LWORK, INFO )
                                 S2 = SECOND( )
                                 TIME = TIME + ( S2-S1 )
                                 IF( INFO.EQ.0 .AND. TIME.LT.TIMMIN )
     $                              GO TO 20
                                 TIMNG( 1 ) = TIME
                                 OPCNT( 1 ) = SASUM( NDATA( 1 ), OPCNT,
     $                                        1 )
                                 CALL SSCAL( NDATA( 1 ), ONE /
     $                                       REAL( NCALL ), OPCNT, 1 )
                                 CALL SSCAL( NDATA( 1 ), ONE /
     $                                       REAL( NCALL ), TIMNG, 1 )
                                 CALL SCOPY( NDATA( 1 ), OPCNT, 1,
     $                                       OPCTBL( 1, ITYPE, NCLS+INB,
     $                                       1 ), 1 )
                                 CALL SCOPY( NDATA( 1 ), TIMNG, 1,
     $                                       TIMTBL( 1, ITYPE, NCLS+INB,
     $                                       1 ), 1 )
                                 DO 30 I = 1, NDATA( 1 )
                                    FLPTBL( I, ITYPE, NCLS+INB, 1 ) = 
     $                              SMFLOP( OPCNT( I ), TIMNG( I ),
     $                                      INFO )
   30                            CONTINUE
   40                         CONTINUE
   50                      CONTINUE
*
                        END IF
*
*                       Generate a matrix of scaling type ISCALE and
*                       rank type IRANK.
*
                        ITYPE = ( IRANK-1 )*3 + ISCALE
                        CALL CQRT15( ISCALE, IRANK, M, N, NRHS, COPYA,
     $                               LDA, COPYB, LDB, COPYS, RANK,
     $                               NORMA, NORMB, ISEED, WORK, LWORK )
*
                        IF( TIMSUB( 2 ) ) THEN
*
*                       Time CGELSX
*
*                       workspace used:
*                       MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)
*
                        LDWORK = MAX( 1, M )
*
*                       CGELSX:  Compute the minimum-norm
*                       solution X to min( norm( A * X - B ) )
*                       using a complete orthogonal factorization.
*
                        NCALL = 0
                        TIME = ZERO
                        CALL SLASET( 'Full', NDATA( 2 ), 1, ZERO, ZERO,
     $                               OPCNT, NDATA( 2 ) )
                        CALL SLASET( 'Full', NDATA( 2 ), 1, ZERO, ZERO,
     $                               TIMNG, NDATA( 2 ) )
   60                   CONTINUE
                        CALL CLACPY( 'Full', M, N, COPYA, LDA,
     $                               A, LDA )
                        CALL CLACPY( 'Full', M, NRHS, COPYB, LDB,
     $                               B, LDB )
                        SRNAMT = 'CGELSX'
                        NCALL = NCALL + 1
                        S1 = SECOND( )
                        CALL CGELSX( M, N, NRHS, A, LDA, B, LDB, IWORK,
     $                               RCOND, CRANK, WORK, RWORK, INFO )
                        S2 = SECOND( )
                        TIME = TIME + ( S2-S1 )
                        IF( INFO.EQ.0 .AND. TIME.LT.TIMMIN )
     $                     GO TO 60
                        TIMNG( 1 ) = TIME
                        OPCNT( 1 ) = SASUM( NDATA( 2 ), OPCNT, 1 )
                        CALL SSCAL( NDATA( 2 ), ONE / REAL( NCALL ),
     $                              OPCNT, 1 )
                        CALL SSCAL( NDATA( 2 ), ONE / REAL( NCALL ),
     $                              TIMNG, 1 )
                        CALL SCOPY( NDATA( 2 ), OPCNT, 1, OPCTBL( 1,
     $                              ITYPE, NCLSX+1, 2 ), 1 )
                        CALL SCOPY( NDATA( 2 ), TIMNG, 1, TIMTBL( 1,
     $                              ITYPE, NCLSX+1, 2 ), 1 )
                        DO 70 I = 1, NDATA( 2 )
                           FLPTBL( I, ITYPE, NCLSX+1, 2 ) = 
     $                     SMFLOP( OPCNT( I ), TIMNG( I ), INFO )
   70                   CONTINUE
*
                        END IF
*
*                       Loop for timing different block sizes.
*
                        DO 140 INB = 1, NNB
                           NB = NBVAL( INB )
                           CALL XLAENV( 1, NB )
                           CALL XLAENV( 3, NXVAL( INB ) )
*
                           IF( TIMSUB( 3 ) ) THEN
*
*                          Time CGELSY
*
*                          CGELSY:  Compute the minimum-norm solution X
*                          to min( norm( A * X - B ) ) using the
*                          rank-revealing orthogonal factorization.
*
*                          Set LWLSY to the adequate value.
*
                           LWLSY = MAX( 1, MNMIN+2*N+NB*( N+1 ),
     $                             2*MNMIN+NB*NRHS )
*
                           NCALL = 0
                           TIME = ZERO
                           CALL SLASET( 'Full', NDATA( 3 ), 1, ZERO,
     $                                  ZERO, OPCNT, NDATA( 3 ) )
                           CALL SLASET( 'Full', NDATA( 3 ), 1, ZERO,
     $                                  ZERO, TIMNG, NDATA( 3 ) )
   80                      CONTINUE
                           CALL CLACPY( 'Full', M, N, COPYA, LDA,
     $                                  A, LDA )
                           CALL CLACPY( 'Full', M, NRHS, COPYB, LDB,
     $                                  B, LDB )
                           SRNAMT = 'CGELSY'
                           NCALL = NCALL + 1
                           S1 = SECOND( )
                           CALL CGELSY( M, N, NRHS, A, LDA, B, LDB,
     $                               IWORK, RCOND, CRANK, WORK, LWLSY,
     $                               RWORK, INFO )
                           S2 = SECOND( )
                           TIME = TIME + ( S2-S1 )
                           IF( INFO.EQ.0 .AND. TIME.LT.TIMMIN )
     $                        GO TO 80
                           TIMNG( 1 ) = TIME
                           OPCNT( 1 ) = SASUM( NDATA( 3 ), OPCNT, 1 )
                           CALL SSCAL( NDATA( 3 ), ONE / REAL( NCALL ),
     $                                 OPCNT, 1 )
                           CALL SSCAL( NDATA( 3 ), ONE / REAL( NCALL ),
     $                                 TIMNG, 1 )
                           CALL SCOPY( NDATA( 3 ), OPCNT, 1, OPCTBL( 1,
     $                                 ITYPE, NCLSY+INB, 3 ), 1 )
                           CALL SCOPY( NDATA( 3 ), TIMNG, 1, TIMTBL( 1,
     $                                 ITYPE, NCLSY+INB, 3 ), 1 )
                           DO 90 I = 1, NDATA( 3 )
                              FLPTBL( I, ITYPE, NCLSY+INB, 3 ) = 
     $                        SMFLOP( OPCNT( I ), TIMNG( I ), INFO )
   90                      CONTINUE
*
                           END IF
*
                           IF( TIMSUB( 4 ) ) THEN
*
*                          Time CGELSS
*
*                          CGELSS:  Compute the minimum-norm solution X
*                          to min( norm( A * X - B ) ) using the SVD.
*
                           NCALL = 0
                           TIME = ZERO
                           CALL SLASET( 'Full', NDATA( 4 ), 1, ZERO,
     $                                  ZERO, OPCNT, NDATA( 4 ) )
                           CALL SLASET( 'Full', NDATA( 4 ), 1, ZERO,
     $                                  ZERO, TIMNG, NDATA( 4 ) )
  100                      CONTINUE
                           CALL CLACPY( 'Full', M, N, COPYA, LDA,
     $                                  A, LDA )
                           CALL CLACPY( 'Full', M, NRHS, COPYB, LDB,
     $                                  B, LDB )
                           SRNAMT = 'CGELSS'
                           NCALL = NCALL + 1
                           S1 = SECOND( )
                           CALL CGELSS( M, N, NRHS, A, LDA, B, LDB,
     $                                  S, RCOND, CRANK, WORK, LWORK,
     $                                  RWORK, INFO )
                           S2 = SECOND( )
                           TIME = TIME + ( S2-S1 )
                           IF( INFO.EQ.0 .AND. TIME.LT.TIMMIN )
     $                        GO TO 100
                           TIMNG( 1 ) = TIME
                           OPCNT( 1 ) = SASUM( NDATA( 4 ), OPCNT, 1 )
                           CALL SSCAL( NDATA( 4 ), ONE / REAL( NCALL ),
     $                                 OPCNT, 1 )
                           CALL SSCAL( NDATA( 4 ), ONE / REAL( NCALL ),
     $                                 TIMNG, 1 )
                           CALL SCOPY( NDATA( 4 ), OPCNT, 1, OPCTBL( 1,
     $                                 ITYPE, NCLSS+INB, 4 ), 1 )
                           CALL SCOPY( NDATA( 4 ), TIMNG, 1, TIMTBL( 1,
     $                                 ITYPE, NCLSS+INB, 4 ), 1 )
                           DO 110 I = 1, NDATA( 4 )
                              FLPTBL( I, ITYPE, NCLSS+INB, 4 ) = 
     $                        SMFLOP( OPCNT( I ), TIMNG( I ), INFO )
  110                      CONTINUE
*
                           END IF
*
                           IF( TIMSUB( 5 ) ) THEN
*
*                          Time CGELSD
*
*                          CGELSD:  Compute the minimum-norm solution X
*                          to min( norm( A * X - B ) ) using a
*                          divide-and-conquer SVD.
*
                           CALL XLAENV( 9, 25 )
                           NCALL = 0
                           TIME = ZERO
                           CALL SLASET( 'Full', NDATA( 5 ), 1, ZERO,
     $                                  ZERO, OPCNT, NDATA( 5 ) )
                           CALL SLASET( 'Full', NDATA( 5 ), 1, ZERO,
     $                                  ZERO, TIMNG, NDATA( 5 ) )
  120                      CONTINUE
                           CALL CLACPY( 'Full', M, N, COPYA, LDA,
     $                                  A, LDA )
                           CALL CLACPY( 'Full', M, NRHS, COPYB, LDB,
     $                                  B, LDB )
                           SRNAMT = 'CGELSD'
                           NCALL = NCALL + 1
                           S1 = SECOND( )
                           CALL CGELSD( M, N, NRHS, A, LDA, B, LDB, S,
     $                                  RCOND, CRANK, WORK, LWORK,
     $                                  RWORK, IWORK, INFO )
                           S2 = SECOND( )
                           TIME = TIME + ( S2-S1 )
                           IF( INFO.EQ.0 .AND. TIME.LT.TIMMIN )
     $                        GO TO 120
                           TIMNG( 1 ) = TIME
                           OPCNT( 1 ) = SASUM( NDATA( 5 ), OPCNT, 1 )
                           CALL SSCAL( NDATA( 5 ), ONE / REAL( NCALL ),
     $                                 OPCNT, 1 )
                           CALL SSCAL( NDATA( 5 ), ONE / REAL( NCALL ),
     $                                 TIMNG, 1 )
                           CALL SCOPY( NDATA( 5 ), OPCNT, 1, OPCTBL( 1,
     $                                 ITYPE, NCLSD+INB, 5 ), 1 )
                           CALL SCOPY( NDATA( 5 ), TIMNG, 1, TIMTBL( 1,
     $                                 ITYPE, NCLSD+INB, 5 ), 1 )
                           DO 130 I = 1, NDATA( 5 )
                              FLPTBL( I, ITYPE, NCLSD+INB, 5 ) = 
     $                        SMFLOP( OPCNT( I ), TIMNG( I ), INFO )
  130                      CONTINUE
*
                           END IF
*
  140                   CONTINUE
  150                CONTINUE
  160             CONTINUE
                  NCLS = NCLS + NNB
                  NCLSY = NCLSY + NNB
                  NCLSS = NCLSS + NNB
                  NCLSD = NCLSD + NNB
  170             CONTINUE
               NCLSX = NCLSX + 1
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
*
*     Print a summary of the results.
*
      DO 220 ISUB = 1, NSUBS
         IF( TIMSUB( ISUB ) ) THEN
            WRITE( NOUT, FMT = 9999 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
            IF( ISUB.EQ.1 ) THEN
               WRITE( NOUT, FMT = 9998 )
            ELSE IF( ISUB.EQ.2 ) THEN
               WRITE( NOUT, FMT = 9997 )
            ELSE IF( ISUB.EQ.3 ) THEN
               WRITE( NOUT, FMT = 9996 )
            ELSE IF( ISUB.EQ.4 ) THEN
               WRITE( NOUT, FMT = 9995 )
            ELSE IF( ISUB.EQ.5 ) THEN
               WRITE( NOUT, FMT = 9994 )
            END IF
            DO 210 ITBL = 1, 3
               IF( ITBL.EQ.1 ) THEN
                  WRITE( NOUT, FMT = 9993 )
                  CALL SPRTLS( ISUB, SUBNAM( ISUB ), NDATA( ISUB ),
     $                         NM, MVAL, NN, NVAL, NNS, NSVAL, NNB,
     $                         NBVAL, NXVAL, NLDA, LDAVAL, MTYPE,
     $                         TIMTBL( 1, 1, 1, ISUB ), NOUT )
               ELSE IF( ITBL.EQ.2 ) THEN
                  WRITE( NOUT, FMT = 9992 )
                  CALL SPRTLS( ISUB, SUBNAM( ISUB ), NDATA( ISUB ),
     $                         NM, MVAL, NN, NVAL, NNS, NSVAL, NNB,
     $                         NBVAL, NXVAL, NLDA, LDAVAL, MTYPE,
     $                         OPCTBL( 1, 1, 1, ISUB ), NOUT )
               ELSE IF( ITBL.EQ.3 ) THEN
                  WRITE( NOUT, FMT = 9991 )
                  CALL SPRTLS( ISUB, SUBNAM( ISUB ), NDATA( ISUB ),
     $                         NM, MVAL, NN, NVAL, NNS, NSVAL, NNB,
     $                         NBVAL, NXVAL, NLDA, LDAVAL, MTYPE,
     $                         FLPTBL( 1, 1, 1, ISUB ), NOUT )
               END IF
  210       CONTINUE
         END IF
  220 CONTINUE
*
  230 CONTINUE
 9999 FORMAT( / / / ' ****** Results for ', A, ' ******' )
 9998 FORMAT( / ' CGELS   : overall performance',
     $        / ' comp. 1 : if M>=N, CGEQRF, QR factorization',
     $        / '           if M< N, CGELQF, QR factorization',
     $        / ' comp. 2 : if M>=N, CUNMQR, multiplication by',
     $          ' reflectors',
     $        / '           if M< N, CUNMLQ, multiplication by',
     $          ' reflectors',
     $        / ' comp. 3 : CTRSM, solution of the triangular',
     $          ' system' /, 
     $        / ' Types 4 to 6 are the conjugate transpose',
     $          ' of types 1 to 3' )
 9997 FORMAT( / ' CGELSX  : overall performance',
     $        / ' comp. 1 : CGEQPF, QR factorization with column',
     $          ' pivoting',
     $        / ' comp. 2 : if RANK<N, CTZRQF, reduction to',
     $          ' triangular form',
     $        / ' comp. 3 : CUNM2R, multiplication by reflectors',
     $        / ' comp. 4 : CTRSM, solution of the triangular',
     $          ' system',
     $        / ' comp. 5 : if RANK<N, CLATZM, multiplication by',
     $          ' reflectors' )
 9996 FORMAT( / ' CGELSY  : overall performance',
     $        / ' comp. 1 : CGEQP3, QR factorization with column',
     $          ' pivoting',
     $        / ' comp. 2 : if RANK<N, CTZRZF, reduction to',
     $          ' triangular form',
     $        / ' comp. 3 : CUNMQR, multiplication by reflectors',
     $        / ' comp. 4 : CTRSM, solution of the triangular',
     $          ' system',
     $        / ' comp. 5 : if RANK<N, CUNMRZ, multiplication by',
     $          ' reflectors' )
 9995 FORMAT( / ' CGELSS: overall performance',
     $        / ' comp. 1 : if M>>N, CGEQRF, QR factorization',
     $        / '                    CUNMQR, multiplication by',
     $          ' reflectors',
     $        / '           if N>>M, CGELQF, QL factorization',
     $        / ' comp. 2 : CGEBRD, reduction to bidiagonal form',
     $        / ' comp. 3 : CUNMBR, multiplication by left',
     $          ' bidiagonalizing vectors',
     $        / '           CUNGBR, generation of right',
     $          ' bidiagonalizing vectors',
     $        / ' comp. 4 : CBDSQR, singular value decomposition',
     $          ' of the bidiagonal matrix',
     $        / ' comp. 5 : multiplication by right bidiagonalizing',
     $          ' vectors',
     $        / '           (CGEMM or CGEMV, and CUNMLQ if N>>M)' )
 9994 FORMAT( / ' CGELSD: overall performance',
     $        / ' comp. 1 : if M>>N, CGEQRF, QR factorization',
     $        / '                    CUNMQR, multiplication by',
     $          ' reflectors',
     $        / '           if N>>M, CGELQF, QL factorization',
     $        / ' comp. 2 : CGEBRD, reduction to bidiagonal form',
     $        / ' comp. 3 : CUNMBR, multiplication by left ',
     $          ' bidiagonalizing vectors',
     $        / '                   multiplication by right',
     $          ' bidiagonalizing vectors',
     $        / ' comp. 4 : CLALSD, singular value decomposition',
     $          ' of the bidiagonal matrix' )
 9993 FORMAT( / / ' *** Time in seconds *** ' )
 9992 FORMAT( / / ' *** Number of floating-point operations *** ' )
 9991 FORMAT( / / ' *** Speed in megaflops *** ' )
      RETURN
*
*     End of CTIMLS
*
      END
