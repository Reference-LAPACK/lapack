      SUBROUTINE DTIMTD( LINE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL,
     $                   NLDA, LDAVAL, TIMMIN, A, B, D, TAU, WORK,
     $                   RESLTS, LDR1, LDR2, LDR3, NOUT )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            LDR1, LDR2, LDR3, NLDA, NM, NN, NNB, NOUT
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      INTEGER            LDAVAL( * ), MVAL( * ), NBVAL( * ), NVAL( * ),
     $                   NXVAL( * )
      DOUBLE PRECISION   A( * ), B( * ), D( * ),
     $                   RESLTS( LDR1, LDR2, LDR3, * ), TAU( * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DTIMTD times the LAPACK routines DSYTRD, DORGTR, and DORMTR and the
*  EISPACK routine TRED1.
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
*  NN      (input) INTEGER
*          The number of values of N contained in the vector NVAL.
*
*  NVAL    (input) INTEGER array, dimension (NN)
*          The values of the matrix column dimension N.
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
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  A       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*          where LDAMAX and NMAX are the maximum values of LDA and N.
*
*  B       (workspace) DOUBLE PRECISION array, dimension (LDAMAX*NMAX)
*
*  D       (workspace) DOUBLE PRECISION array, dimension (2*NMAX-1)
*
*  TAU     (workspace) DOUBLE PRECISION array, dimension (NMAX)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (NMAX*NBMAX)
*          where NBMAX is the maximum value of NB.
*
*  RESLTS  (workspace) DOUBLE PRECISION array, dimension
*                      (LDR1,LDR2,LDR3,4*NN+3)
*          The timing results for each subroutine over the relevant
*          values of M, (NB,NX), LDA, and N.
*
*  LDR1    (input) INTEGER
*          The first dimension of RESLTS.  LDR1 >= max(1,NNB).
*
*  LDR2    (input) INTEGER
*          The second dimension of RESLTS.  LDR2 >= max(1,NM).
*
*  LDR3    (input) INTEGER
*          The third dimension of RESLTS.  LDR3 >= max(1,2*NLDA).
*
*  NOUT    (input) INTEGER
*          The unit number for output.
*
*  Internal Parameters
*  ===================
*
*  MODE    INTEGER
*          The matrix type.  MODE = 3 is a geometric distribution of
*          eigenvalues.  See DLATMS for further details.
*
*  COND    DOUBLE PRECISION
*          The condition number of the matrix.  The singular values are
*          set to values from DMAX to DMAX/COND.
*
*  DMAX    DOUBLE PRECISION
*          The magnitude of the largest singular value.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NSUBS
      PARAMETER          ( NSUBS = 4 )
      INTEGER            MODE
      DOUBLE PRECISION   COND, DMAX
      PARAMETER          ( MODE = 3, COND = 100.0D0, DMAX = 1.0D0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          LAB1, LAB2, SIDE, TRANS, UPLO
      CHARACTER*3        PATH
      CHARACTER(32)      CNAME
      INTEGER            I, I3, I4, IC, ICL, ILDA, IM, IN, INB, INFO,
     $                   ISIDE, ISUB, ITOFF, ITRAN, IUPLO, LDA, LW, M,
     $                   M1, N, N1, NB, NX
      DOUBLE PRECISION   OPS, S1, S2, TIME, UNTIME
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER          SIDES( 2 ), TRANSS( 2 ), UPLOS( 2 )
      CHARACTER(32)      SUBNAM( NSUBS )
      INTEGER            ISEED( 4 ), RESEED( 4 )
*     ..
*     .. External Functions ..
      INTEGER ILA_LEN_TRIM
      EXTERNAL ILA_LEN_TRIM
      DOUBLE PRECISION   DMFLOP, DOPLA, DSECND
      EXTERNAL           DMFLOP, DOPLA, DSECND
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMCK, ATIMIN, DLACPY, DLATMS, DORGTR, DORMTR,
     $                   DPRTB3, DPRTBL, DSYTRD, DTIMMG, ICOPY, TRED1,
     $                   XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'DSYTRD', 'TRED1', 'DORGTR',
     $                   'DORMTR' /
      DATA               SIDES / 'L', 'R' / , TRANSS / 'N', 'T' / ,
     $                   UPLOS / 'U', 'L' /
      DATA               ISEED / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
*     Extract the timing request from the input line.
*
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'TD'
      CALL ATIMIN( PATH, LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
      IF( INFO.NE.0 )
     $   GO TO 240
*
*     Check that M <= LDA for the input values.
*
      CNAME = LINE( 1: 6 )
      CALL ATIMCK( 2, CNAME, NM, MVAL, NLDA, LDAVAL, NOUT, INFO )
      IF( INFO.GT.0 ) THEN
         WRITE( NOUT, FMT = 9999 )CNAME(1:ILA_LEN_TRIM(CNAME))
         GO TO 240
      END IF
*
*     Check that K <= LDA for DORMTR
*
      IF( TIMSUB( 4 ) ) THEN
         CALL ATIMCK( 3, CNAME, NN, NVAL, NLDA, LDAVAL, NOUT, INFO )
         IF( INFO.GT.0 ) THEN
            WRITE( NOUT, FMT = 9999 )
     $     SUBNAM( 4 )(1:ILA_LEN_TRIM( SUBNAM( 4 ) ))
            TIMSUB( 4 ) = .FALSE.
         END IF
      END IF
*
*     Do first for UPLO = 'U', then for UPLO = 'L'
*
      DO 150 IUPLO = 1, 2
         UPLO = UPLOS( IUPLO )
*
*        Do for each value of M:
*
         DO 140 IM = 1, NM
            M = MVAL( IM )
            CALL ICOPY( 4, ISEED, 1, RESEED, 1 )
*
*           Do for each value of LDA:
*
            DO 130 ILDA = 1, NLDA
               LDA = LDAVAL( ILDA )
               I3 = ( IUPLO-1 )*NLDA + ILDA
*
*              Do for each pair of values (NB, NX) in NBVAL and NXVAL.
*
               DO 120 INB = 1, NNB
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
                  NX = NXVAL( INB )
                  CALL XLAENV( 3, NX )
                  LW = MAX( 1, M*MAX( 1, NB ) )
*
*                 Generate a test matrix of order M.
*
                  CALL ICOPY( 4, RESEED, 1, ISEED, 1 )
                  CALL DLATMS( M, M, 'Uniform', ISEED, 'Symmetric', TAU,
     $                         MODE, COND, DMAX, M, M, 'No packing', B,
     $                         LDA, WORK, INFO )
*
                  IF( TIMSUB( 2 ) .AND. INB.EQ.1 .AND. IUPLO.EQ.2 ) THEN
*
*                    TRED1:  Eispack reduction using orthogonal
*                    transformations.
*
                     CALL DLACPY( UPLO, M, M, B, LDA, A, LDA )
                     IC = 0
                     S1 = DSECND( )
   10                CONTINUE
                     CALL TRED1( LDA, M, A, D, D( M+1 ), D( M+1 ) )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL DLACPY( UPLO, M, M, B, LDA, A, LDA )
                        GO TO 10
                     END IF
*
*                    Subtract the time used in DLACPY.
*
                     ICL = 1
                     S1 = DSECND( )
   20                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL DLACPY( UPLO, M, M, B, LDA, A, LDA )
                        GO TO 20
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPLA( 'DSYTRD', M, M, -1, -1, NB )
                     RESLTS( INB, IM, ILDA, 2 ) = DMFLOP( OPS, TIME,
     $                  INFO )
                  END IF
*
                  IF( TIMSUB( 1 ) ) THEN
*
*                    DSYTRD:  Reduction to tridiagonal form
*
                     CALL DLACPY( UPLO, M, M, B, LDA, A, LDA )
                     IC = 0
                     S1 = DSECND( )
   30                CONTINUE
                     CALL DSYTRD( UPLO, M, A, LDA, D, D( M+1 ), TAU,
     $                            WORK, LW, INFO )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL DLACPY( UPLO, M, M, B, LDA, A, LDA )
                        GO TO 30
                     END IF
*
*                    Subtract the time used in DLACPY.
*
                     ICL = 1
                     S1 = DSECND( )
   40                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL DLACPY( UPLO, M, M, A, LDA, B, LDA )
                        GO TO 40
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
                     OPS = DOPLA( 'DSYTRD', M, M, -1, -1, NB )
                     RESLTS( INB, IM, I3, 1 ) = DMFLOP( OPS, TIME,
     $                  INFO )
                  ELSE
*
*                    If DSYTRD was not timed, generate a matrix and
*                    factor it using DSYTRD anyway so that the factored
*                    form of the matrix can be used in timing the other
*                    routines.
*
                     CALL DLACPY( UPLO, M, M, B, LDA, A, LDA )
                     CALL DSYTRD( UPLO, M, A, LDA, D, D( M+1 ), TAU,
     $                            WORK, LW, INFO )
                  END IF
*
                  IF( TIMSUB( 3 ) ) THEN
*
*                    DORGTR:  Generate the orthogonal matrix Q from the
*                    reduction to Hessenberg form A = Q*H*Q'
*
                     CALL DLACPY( UPLO, M, M, A, LDA, B, LDA )
                     IC = 0
                     S1 = DSECND( )
   50                CONTINUE
                     CALL DORGTR( UPLO, M, B, LDA, TAU, WORK, LW, INFO )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN ) THEN
                        CALL DLACPY( UPLO, M, M, A, LDA, B, LDA )
                        GO TO 50
                     END IF
*
*                    Subtract the time used in DLACPY.
*
                     ICL = 1
                     S1 = DSECND( )
   60                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     ICL = ICL + 1
                     IF( ICL.LE.IC ) THEN
                        CALL DLACPY( UPLO, M, M, A, LDA, B, LDA )
                        GO TO 60
                     END IF
*
                     TIME = ( TIME-UNTIME ) / DBLE( IC )
*
*                    Op count for DORGTR:  same as
*                       DORGQR( N-1, N-1, N-1, ... )
*
                     OPS = DOPLA( 'DORGQR', M-1, M-1, M-1, -1, NB )
                     RESLTS( INB, IM, I3, 3 ) = DMFLOP( OPS, TIME,
     $                  INFO )
                  END IF
*
                  IF( TIMSUB( 4 ) ) THEN
*
*                    DORMTR:  Multiply by Q stored as a product of
*                    elementary transformations
*
                     I4 = 3
                     DO 110 ISIDE = 1, 2
                        SIDE = SIDES( ISIDE )
                        DO 100 IN = 1, NN
                           N = NVAL( IN )
                           LW = MAX( 1, MAX( 1, NB )*N )
                           IF( ISIDE.EQ.1 ) THEN
                              M1 = M
                              N1 = N
                           ELSE
                              M1 = N
                              N1 = M
                           END IF
                           ITOFF = 0
                           DO 90 ITRAN = 1, 2
                              TRANS = TRANSS( ITRAN )
                              CALL DTIMMG( 0, M1, N1, B, LDA, 0, 0 )
                              IC = 0
                              S1 = DSECND( )
   70                         CONTINUE
                              CALL DORMTR( SIDE, UPLO, TRANS, M1, N1, A,
     $                                     LDA, TAU, B, LDA, WORK, LW,
     $                                     INFO )
                              S2 = DSECND( )
                              TIME = S2 - S1
                              IC = IC + 1
                              IF( TIME.LT.TIMMIN ) THEN
                                 CALL DTIMMG( 0, M1, N1, B, LDA, 0, 0 )
                                 GO TO 70
                              END IF
*
*                             Subtract the time used in DTIMMG.
*
                              ICL = 1
                              S1 = DSECND( )
   80                         CONTINUE
                              S2 = DSECND( )
                              UNTIME = S2 - S1
                              ICL = ICL + 1
                              IF( ICL.LE.IC ) THEN
                                 CALL DTIMMG( 0, M1, N1, B, LDA, 0, 0 )
                                 GO TO 80
                              END IF
*
                              TIME = ( TIME-UNTIME ) / DBLE( IC )
*
*                             Op count for DORMTR, SIDE='L':  same as
*                                DORMQR( 'L', TRANS, M-1, N, M-1, ...)
*
*                             Op count for DORMTR, SIDE='R':  same as
*                                DORMQR( 'R', TRANS, M, N-1, N-1, ...)
*
                              IF( ISIDE.EQ.1 ) THEN
                                 OPS = DOPLA( 'DORMQR', M1-1, N1, M1-1,
     $                                 -1, NB )
                              ELSE
                                 OPS = DOPLA( 'DORMQR', M1, N1-1, N1-1,
     $                                 1, NB )
                              END IF
*
                              RESLTS( INB, IM, I3,
     $                           I4+ITOFF+IN ) = DMFLOP( OPS, TIME,
     $                           INFO )
                              ITOFF = NN
   90                      CONTINUE
  100                   CONTINUE
                        I4 = I4 + 2*NN
  110                CONTINUE
                  END IF
*
  120          CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
*
*     Print tables of results for DSYTRD, TRED1, and DORGTR
*
      DO 180 ISUB = 1, NSUBS - 1
         IF( .NOT.TIMSUB( ISUB ) )
     $      GO TO 180
         WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
         IF( NLDA.GT.1 ) THEN
            DO 160 I = 1, NLDA
               WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  160       CONTINUE
         END IF
         IF( ISUB.EQ.2 ) THEN
            WRITE( NOUT, FMT = * )
            CALL DPRTB3( ' ', 'N', 1, NBVAL, NXVAL, NM, MVAL, NLDA,
     $                   RESLTS( 1, 1, 1, ISUB ), LDR1, LDR2, NOUT )
         ELSE
            I3 = 1
            DO 170 IUPLO = 1, 2
               WRITE( NOUT, FMT = 9996 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )),
     $     UPLOS( IUPLO )
               CALL DPRTB3( '(  NB,  NX)', 'N', NNB, NBVAL, NXVAL, NM,
     $                      MVAL, NLDA, RESLTS( 1, 1, I3, ISUB ), LDR1,
     $                      LDR2, NOUT )
               I3 = I3 + NLDA
  170       CONTINUE
         END IF
  180 CONTINUE
*
*     Print tables of results for DORMTR
*
      ISUB = 4
      IF( TIMSUB( ISUB ) ) THEN
         I4 = 3
         DO 230 ISIDE = 1, 2
            IF( ISIDE.EQ.1 ) THEN
               LAB1 = 'M'
               LAB2 = 'N'
               IF( NLDA.GT.1 ) THEN
                  WRITE( NOUT, FMT = 9998 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) ))
                  DO 190 I = 1, NLDA
                     WRITE( NOUT, FMT = 9997 )I, LDAVAL( I )
  190             CONTINUE
               END IF
            ELSE
               LAB1 = 'N'
               LAB2 = 'M'
            END IF
            DO 220 ITRAN = 1, 2
               DO 210 IN = 1, NN
                  I3 = 1
                  DO 200 IUPLO = 1, 2
                     WRITE( NOUT, FMT = 9995 )
     $     SUBNAM( ISUB )(1:ILA_LEN_TRIM( SUBNAM( ISUB ) )),
     $                  SIDES( ISIDE ), UPLOS( IUPLO ), TRANSS( ITRAN ),
     $                  LAB2, NVAL( IN )
                     CALL DPRTBL( 'NB', LAB1, NNB, NBVAL, NM, MVAL,
     $                            NLDA, RESLTS( 1, 1, I3, I4+IN ), LDR1,
     $                            LDR2, NOUT )
                     I3 = I3 + NLDA
  200             CONTINUE
  210          CONTINUE
               I4 = I4 + NN
  220       CONTINUE
  230    CONTINUE
      END IF
  240 CONTINUE
*
*     Print a table of results for each timed routine.
*
 9999 FORMAT( 1X, A, ' timing run not attempted', / )
 9998 FORMAT( / ' *** Speed of ', A, ' in megaflops *** ' )
 9997 FORMAT( 5X, 'line ', I2, ' with LDA = ', I5 )
 9996 FORMAT( / 5X, A, ' with UPLO = ''', A1, '''', / )
 9995 FORMAT( / 5X, A, ' with SIDE = ''', A1, ''', UPLO = ''', A1,
     $      ''', TRANS = ''', A1, ''', ', A1, ' =', I6, / )
      RETURN
*
*     End of DTIMTD
*
      END
