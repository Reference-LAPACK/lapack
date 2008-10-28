      PROGRAM CTIMAA
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*  Purpose
*  =======
*
*  CTIMAA is the timing program for the COMPLEX LAPACK
*  routines.  This program collects performance data for the factor,
*  solve, and inverse routines used in solving systems of linear
*  equations, and also for the orthogonal factorization and reduction
*  routines used in solving least squares problems and matrix eigenvalue
*  problems.
*
*  The subprograms call a REAL function SECOND with no
*  arguments which is assumed to return the central-processor time in
*  seconds from some fixed starting time.
*
*  The program is driven by a short data file, which specifies values
*  for the matrix dimensions M, N and K, for the blocking parameters
*  NB and NX, and for the leading array dimension LDA.  A minimum time
*  for each subroutine is included for timing small problems or for
*  obtaining results on a machine with an inaccurate SECOND function.
*
*  The matrix dimensions M, N, and K correspond to the three dimensions
*  m, n, and k in the Level 3 BLAS.  When timing the LAPACK routines for
*  square matrices, M and N correspond to the matrix dimensions m and n,
*  and K is the number of right-hand sides (nrhs) for the solves.  When
*  timing the LAPACK routines for band matrices, M is the matrix order
*  m, N is the half-bandwidth (kl, ku, or kd in the LAPACK notation),
*  and K is again the number of right-hand sides.
*
*  The first 13 records of the data file are read using list-directed
*  input.  The first line of input is printed as the first line of
*  output and can be used to identify different sets of results.  To
*  assist with debugging an input file, the values are printed out as
*  they are read in.
*
*  The following records are read using the format (A).  For these
*  records, the first 6 characters are reserved for the path or
*  subroutine name.  If a path name is used, the characters after the
*  path name indicate the routines in the path to be timed, where
*  'T' or 't' means 'Time this routine'.  If the line is blank after the
*  path name, all routines in the path are timed.  If fewer characters
*  appear than routines in a path, the remaining characters are assumed
*  to be 'F'.  For example, the following 3 lines are equivalent ways of
*  requesting timing of CGETRF:
*  CGE    T F F
*  CGE    T
*  CGETRF
*
*  An annotated example of a data file can be obtained by deleting the
*  first 3 characters from the following 32 lines:
*  LAPACK timing, COMPLEX square matrices
*  5                                Number of values of M
*  100 200 300 400 500              Values of M (row dimension)
*  5                                Number of values of N
*  100 200 300 400 500              Values of N (column dimension)
*  2                                Number of values of K
*  100 400                          Values of K
*  5                                Number of values of NB
*  1 16  32  48  64                 Values of NB (blocksize)
*  0 48 128 128 128                 Values of NX (crossover point)
*  2                                Number of values of LDA
*  512 513                          Values of LDA (leading dimension)
*  0.0                              Minimum time in seconds
*  CGE    T T T
*  CPO    T T T
*  CPP    T T T
*  CHE    T T T
*  CHP    T T T
*  CSY    T T T
*  CSP    T T T
*  CTR    T T
*  CTP    T T
*  CQR    T T F
*  CLQ    T T F
*  CQL    T T F
*  CRQ    T T F
*  CQP    T
*  CHR    T T F F
*  CTD    T T F F
*  CBR    T F F
*  CLS    T T T T T T
*
*  The routines are timed for all combinations of applicable values of
*  M, N, K, NB, NX, and LDA, and for all combinations of options such as
*  UPLO and TRANS.  For Level 2 BLAS timings, values of NB are used for
*  INCX.  Certain subroutines, such as the QR factorization, treat the
*  values of M and N as ordered pairs and operate on M x N matrices.
*
*  Internal Parameters
*  ===================
*
*  NMAX    INTEGER
*          The maximum value of M or N for square matrices.
*
*  LDAMAX  INTEGER
*          The maximum value of LDA.
*
*  NMAXB   INTEGER
*          The maximum value of N for band matrices.
*
*  MAXVAL  INTEGER
*          The maximum number of values that can be read in for M, N,
*          K, NB, or NX.
*
*  MXNLDA  INTEGER
*          The maximum number of values that can be read in for LDA.
*
*  NIN     INTEGER
*          The unit number for input.  Currently set to 5 (std input).
*
*  NOUT    INTEGER
*          The unit number for output.  Currently set to 6 (std output).
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NMAX, LDAMAX, NMAXB
      PARAMETER          ( NMAX = 512, LDAMAX = NMAX+20, NMAXB = 5000 )
      INTEGER            LA
      PARAMETER          ( LA = NMAX*LDAMAX )
      INTEGER            MAXVAL, MXNLDA
      PARAMETER          ( MAXVAL = 12, MXNLDA = 4 )
      INTEGER            MAXPRM
      PARAMETER          ( MAXPRM = MXNLDA*(MAXVAL+1) )
      INTEGER            MAXSZS
      PARAMETER          ( MAXSZS = MAXVAL*MAXVAL*MAXVAL )
      INTEGER            NIN, NOUT
      PARAMETER          ( NIN = 5, NOUT = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            BLAS, LDAMOK, LDANOK, LDAOK, MOK, NOK, NXNBOK
      CHARACTER          C1
      CHARACTER*2        C2
      CHARACTER*3        C3
      CHARACTER*80       LINE
      INTEGER            I, I2, J2, L, LDR1, LDR2, LDR3, MAXK, MAXLDA,
     $                   MAXM, MAXN, MAXNB, MKMAX, NEED, NK, NLDA, NM,
     $                   NN, NNB
      REAL               S1, S2, TIMMIN
*     ..
*     .. Local Arrays ..
      INTEGER            IWORK( 2*NMAXB ), KVAL( MAXVAL ),
     $                   LDAVAL( MXNLDA ), MVAL( MAXVAL ),
     $                   NBVAL( MAXVAL ), NVAL( MAXVAL ),
     $                   NXVAL( MAXVAL )
      REAL               D( 2*NMAX ),
     $                   FLPTBL( 6*6*MAXSZS*MAXPRM*5 ),
     $                   OPCTBL( 6*6*MAXSZS*MAXPRM*5 ),
     $                   RESLTS( MAXVAL, MAXVAL, 2*MXNLDA, 4*MAXVAL ),
     $                   RWORK( 150*NMAX+2*MAXVAL ), S( NMAX*2 ),
     $                   TIMTBL( 6*6*MAXSZS*MAXPRM*5 )
      COMPLEX            A( LA, 3 ), B( LA, 3 ), E( 2*NMAX ),
     $                   WORK( NMAX, NMAX+MAXVAL+10 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, LSAMEN
      REAL               SECOND
      EXTERNAL           LSAME, LSAMEN, SECOND
*     ..
*     .. External Subroutines ..
      EXTERNAL           CTIMB2, CTIMB3, CTIMBR, CTIMGB, CTIMGE, CTIMGT,
     $                   CTIMHE, CTIMHP, CTIMHR, CTIMLQ, CTIMLS, CTIMMM,
     $                   CTIMMV, CTIMPB, CTIMPO, CTIMPP, CTIMPT, CTIMQ3,
     $                   CTIMQL, CTIMQP, CTIMQR, CTIMRQ, CTIMSP, CTIMSY,
     $                   CTIMTB, CTIMTD, CTIMTP, CTIMTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Arrays in Common ..
      INTEGER            IPARMS( 100 )
*     ..
*     .. Common blocks ..
      COMMON             / CLAENV / IPARMS
*     ..
*     .. Executable Statements ..
*
      S1 = SECOND( )
      LDR1 = MAXVAL
      LDR2 = MAXVAL
      LDR3 = 2*MXNLDA
      WRITE( NOUT, FMT = 9983 )
*
*     Read the first line.  The first four characters must be 'BLAS'
*     for the BLAS data file format to be used.  Otherwise, the LAPACK
*     data file format is assumed.
*
      READ( NIN, FMT = '( A80 )' )LINE
      BLAS = LSAMEN( 4, LINE, 'BLAS' )
*
*     Find the last non-blank and print the first line of input as the
*     first line of output.
*
      DO 10 L = 80, 1, -1
         IF( LINE( L: L ).NE.' ' )
     $      GO TO 20
   10 CONTINUE
      L = 1
   20 CONTINUE
      WRITE( NOUT, FMT = '( 1X, A, / )' )LINE( 1: L )
      WRITE( NOUT, FMT = 9992 )
*
*     Read in NM and the values for M.
*
      READ( NIN, FMT = * )NM
      IF( NM.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'M', 'NM', MAXVAL
         NM = MAXVAL
      END IF
      READ( NIN, FMT = * )( MVAL( I ), I = 1, NM )
      WRITE( NOUT, FMT = 9991 )'M:     ', ( MVAL( I ), I = 1, NM )
*
*     Check that  M <= NMAXB for all values of M.
*
      MOK = .TRUE.
      MAXM = 0
      DO 30 I = 1, NM
         MAXM = MAX( MVAL( I ), MAXM )
         IF( MVAL( I ).GT.NMAXB ) THEN
            WRITE( NOUT, FMT = 9997 )'M', MVAL( I ), NMAXB
            MOK = .FALSE.
         END IF
   30 CONTINUE
      IF( .NOT.MOK )
     $   WRITE( NOUT, FMT = * )
*
*     Read in NN and the values for N.
*
      READ( NIN, FMT = * )NN
      IF( NN.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'N', 'NN', MAXVAL
         NN = MAXVAL
      END IF
      READ( NIN, FMT = * )( NVAL( I ), I = 1, NN )
      WRITE( NOUT, FMT = 9991 )'N:     ', ( NVAL( I ), I = 1, NN )
*
*     Check that  N <= NMAXB for all values of N.
*
      NOK = .TRUE.
      MAXN = 0
      DO 40 I = 1, NN
         MAXN = MAX( NVAL( I ), MAXN )
         IF( NVAL( I ).GT.NMAXB ) THEN
            WRITE( NOUT, FMT = 9997 )'N', NVAL( I ), NMAXB
            NOK = .FALSE.
         END IF
   40 CONTINUE
      IF( .NOT.NOK )
     $   WRITE( NOUT, FMT = * )
*
*     Read in NK and the values for K.
*
      READ( NIN, FMT = * )NK
      IF( NK.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'K', 'NK', MAXVAL
         NK = MAXVAL
      END IF
      READ( NIN, FMT = * )( KVAL( I ), I = 1, NK )
      WRITE( NOUT, FMT = 9991 )'K:     ', ( KVAL( I ), I = 1, NK )
*
*     Find the maximum value of K (= NRHS).
*
      MAXK = 0
      DO 50 I = 1, NK
         MAXK = MAX( KVAL( I ), MAXK )
   50 CONTINUE
      MKMAX = MAXM*MAX( 2, MAXK )
*
*     Read in NNB and the values for NB.  For the BLAS input files,
*     NBVAL is used to store values for INCX and INCY.
*
      READ( NIN, FMT = * )NNB
      IF( NNB.GT.MAXVAL ) THEN
         WRITE( NOUT, FMT = 9999 )'NB', 'NNB', MAXVAL
         NNB = MAXVAL
      END IF
      READ( NIN, FMT = * )( NBVAL( I ), I = 1, NNB )
*
*     Find the maximum value of NB.
*
      MAXNB = 0
      DO 60 I = 1, NNB
         MAXNB = MAX( NBVAL( I ), MAXNB )
   60 CONTINUE
*
      IF( BLAS ) THEN
         WRITE( NOUT, FMT = 9991 )'INCX:  ', ( NBVAL( I ), I = 1, NNB )
         DO 70 I = 1, NNB
            NXVAL( I ) = 0
   70    CONTINUE
      ELSE
*
*        LAPACK data files:  Read in the values for NX.
*
         READ( NIN, FMT = * )( NXVAL( I ), I = 1, NNB )
*
         WRITE( NOUT, FMT = 9991 )'NB:    ', ( NBVAL( I ), I = 1, NNB )
         WRITE( NOUT, FMT = 9991 )'NX:    ', ( NXVAL( I ), I = 1, NNB )
      END IF
*
*     Read in NLDA and the values for LDA.
*
      READ( NIN, FMT = * )NLDA
      IF( NLDA.GT.MXNLDA ) THEN
         WRITE( NOUT, FMT = 9999 )'LDA', 'NLDA', MXNLDA
         NLDA = MXNLDA
      END IF
      READ( NIN, FMT = * )( LDAVAL( I ), I = 1, NLDA )
      WRITE( NOUT, FMT = 9991 )'LDA:   ', ( LDAVAL( I ), I = 1, NLDA )
*
*     Check that LDA >= 1 for all values of LDA.
*
      LDAOK = .TRUE.
      MAXLDA = 0
      DO 80 I = 1, NLDA
         MAXLDA = MAX( LDAVAL( I ), MAXLDA )
         IF( LDAVAL( I ).LE.0 ) THEN
            WRITE( NOUT, FMT = 9998 )LDAVAL( I )
            LDAOK = .FALSE.
         END IF
   80 CONTINUE
      IF( .NOT.LDAOK )
     $   WRITE( NOUT, FMT = * )
*
*     Check that MAXLDA*MAXN <= LA (for the dense routines).
*
      LDANOK = .TRUE.
      NEED = MAXLDA*MAXN
      IF( NEED.GT.LA ) THEN
         WRITE( NOUT, FMT = 9995 )MAXLDA, MAXN, NEED
         LDANOK = .FALSE.
      END IF
*
*     Check that MAXLDA*MAXM + MAXM*MAXK <= 3*LA (for band routines).
*
      LDAMOK = .TRUE.
      NEED = MAXLDA*MAXM + MAXM*MAXK
      IF( NEED.GT.3*LA ) THEN
         NEED = ( NEED+2 ) / 3
         WRITE( NOUT, FMT = 9994 )MAXLDA, MAXM, MAXK, NEED
         LDAMOK = .FALSE.
      END IF
*
*     Check that MAXN*MAXNB (or MAXN*INCX) <= LA.
*
      NXNBOK = .TRUE.
      NEED = MAXN*MAXNB
      IF( NEED.GT.LA ) THEN
         WRITE( NOUT, FMT = 9996 )MAXN, MAXNB, NEED
         NXNBOK = .FALSE.
      END IF
*
      IF( .NOT.( MOK .AND. NOK .AND. LDAOK .AND. LDANOK .AND. NXNBOK ) )
     $     THEN
         WRITE( NOUT, FMT = 9984 )
         GO TO 110
      END IF
      IF( .NOT.LDAMOK )
     $   WRITE( NOUT, FMT = * )
*
*     Read the minimum time to time a subroutine.
*
      WRITE( NOUT, FMT = * )
      READ( NIN, FMT = * )TIMMIN
      WRITE( NOUT, FMT = 9993 )TIMMIN
      WRITE( NOUT, FMT = * )
*
*     Read the first input line.
*
      READ( NIN, FMT = '(A)', END = 100 )LINE
*
*     If the first record is the special signal 'NONE', then get the
*     next line but don't time CGEMV and CGEMM.
*
      IF( LSAMEN( 4, LINE, 'NONE' ) ) THEN
         READ( NIN, FMT = '(A)', END = 100 )LINE
      ELSE
         WRITE( NOUT, FMT = 9990 )
*
*        If the first record is the special signal 'BAND', then time
*        the band routine CGBMV and CGEMM with N = K.
*
         IF( LSAMEN( 4, LINE, 'BAND' ) ) THEN
            IF( LDAMOK ) THEN
               IF( MKMAX.GT.LA ) THEN
                  I2 = 2*LA - MKMAX + 1
                  J2 = 2
               ELSE
                  I2 = LA - MKMAX + 1
                  J2 = 3
               END IF
               CALL CTIMMV( 'CGBMV ', NM, MVAL, NN, NVAL, NLDA, LDAVAL,
     $                      TIMMIN, A( 1, 1 ), MKMAX / 2, A( I2, J2 ),
     $                      A( LA-MKMAX / 2+1, 3 ), RESLTS, LDR1, LDR2,
     $                      NOUT )
            ELSE
               WRITE( NOUT, FMT = 9989 )'CGBMV '
            END IF
            CALL CTIMMM( 'CGEMM ', 'K', NN, NVAL, NLDA, LDAVAL, TIMMIN,
     $                   A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), RESLTS, LDR1,
     $                   LDR2, NOUT )
            READ( NIN, FMT = '(A)', END = 100 )LINE
*
         ELSE
*
*           Otherwise time CGEMV and CGEMM.
*
            CALL CTIMMV( 'CGEMV ', NN, NVAL, NNB, NBVAL, NLDA, LDAVAL,
     $                   TIMMIN, A( 1, 1 ), LA, A( 1, 2 ), A( 1, 3 ),
     $                   RESLTS, LDR1, LDR2, NOUT )
            CALL CTIMMM( 'CGEMM ', 'N', NN, NVAL, NLDA, LDAVAL, TIMMIN,
     $                   A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), RESLTS, LDR1,
     $                   LDR2, NOUT )
         END IF
      END IF
*
*     Call the appropriate timing routine for each input line.
*
      WRITE( NOUT, FMT = 9988 )
   90 CONTINUE
      C1 = LINE( 1: 1 )
      C2 = LINE( 2: 3 )
      C3 = LINE( 4: 6 )
*
*     Check first character for correct precision.
*
      IF( .NOT.LSAME( C1, 'Complex precision' ) ) THEN
         WRITE( NOUT, FMT = 9987 )LINE( 1: 6 )
*
      ELSE IF( LSAMEN( 2, C2, 'B2' ) .OR. LSAMEN( 3, C3, 'MV ' ) .OR.
     $         LSAMEN( 3, C3, 'SV ' ) .OR. LSAMEN( 3, C3, 'R  ' ) .OR.
     $         LSAMEN( 3, C3, 'RC ' ) .OR. LSAMEN( 3, C3, 'RU ' ) .OR.
     $         LSAMEN( 3, C3, 'R2 ' ) ) THEN
*
*        Level 2 BLAS
*
         CALL CTIMB2( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NNB, NBVAL,
     $                NLDA, LDAVAL, LA, TIMMIN, A( 1, 1 ), A( 1, 2 ),
     $                A( 1, 3 ), RESLTS, LDR1, LDR2, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'B3' ) .OR. LSAMEN( 3, C3, 'MM ' ) .OR.
     $         LSAMEN( 3, C3, 'SM ' ) .OR. LSAMEN( 3, C3, 'RK ' ) .OR.
     $         LSAMEN( 3, C3, 'R2K' ) ) THEN
*
*        Level 3 BLAS
*
         CALL CTIMB3( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NLDA, LDAVAL,
     $                TIMMIN, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ), RESLTS,
     $                LDR1, LDR2, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'QR' ) .OR. LSAMEN( 2, C3, 'QR' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'QR' ) ) THEN
*
*        QR routines
*
         CALL CTIMQR( LINE, NN, MVAL, NVAL, NK, KVAL, NNB, NBVAL, NXVAL,
     $                NLDA, LDAVAL, TIMMIN, A( 1, 1 ), E, A( 1, 2 ),
     $                A( 1, 3 ), D, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'LQ' ) .OR. LSAMEN( 2, C3, 'LQ' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'LQ' ) ) THEN
*
*        LQ routines
*
         CALL CTIMLQ( LINE, NN, MVAL, NVAL, NK, KVAL, NNB, NBVAL, NXVAL,
     $                NLDA, LDAVAL, TIMMIN, A( 1, 1 ), E, A( 1, 2 ),
     $                A( 1, 3 ), D, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'QL' ) .OR. LSAMEN( 2, C3, 'QL' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'QL' ) ) THEN
*
*        QL routines
*
         CALL CTIMQL( LINE, NN, MVAL, NVAL, NK, KVAL, NNB, NBVAL, NXVAL,
     $                NLDA, LDAVAL, TIMMIN, A( 1, 1 ), E, A( 1, 2 ),
     $                A( 1, 3 ), D, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'RQ' ) .OR. LSAMEN( 2, C3, 'RQ' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'RQ' ) ) THEN
*
*        RQ routines
*
         CALL CTIMRQ( LINE, NN, MVAL, NVAL, NK, KVAL, NNB, NBVAL, NXVAL,
     $                NLDA, LDAVAL, TIMMIN, A( 1, 1 ), E, A( 1, 2 ),
     $                A( 1, 3 ), D, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'QP' ) .OR. LSAMEN( 3, C3, 'QPF' ) ) THEN
*
*        QR with column pivoting
*
         CALL CTIMQP( LINE, NM, MVAL, NVAL, NLDA, LDAVAL, TIMMIN,
     $                A( 1, 1 ), A( 1, 2 ), E, A( 1, 3 ), D, IWORK,
     $                RESLTS, LDR1, LDR2, NOUT )
*
*        Rank-Revealing QR factorization
*
         CALL CTIMQ3( LINE, NM, MVAL, NVAL, NNB, NBVAL, NXVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), E,
     $                A( 1, 3 ), D, IWORK, RESLTS, LDR1, LDR2, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'HR' ) .OR. LSAMEN( 3, C3, 'HRD' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'HR' ) ) THEN
*
*        Reduction to Hessenberg form
*
         CALL CTIMHR( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NXVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), E, A( 1, 2 ),
     $                A( 1, 3 ), D, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'TD' ) .OR. LSAMEN( 3, C3, 'TRD' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'TR' ) ) THEN
*
*        Reduction to tridiagonal form
*
         CALL CTIMTD( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NXVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), D, E,
     $                A( 1, 3 ), RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'BR' ) .OR. LSAMEN( 3, C3, 'BRD' ) .OR.
     $         LSAMEN( 2, C3( 2: 3 ), 'BR' ) ) THEN
*
*        Reduction to bidiagonal form
*
         CALL CTIMBR( LINE, NN, MVAL, NVAL, NK, KVAL, NNB, NBVAL, NXVAL,
     $                NLDA, LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), D, E,
     $                A( 1, 3 ), RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'GE' ) ) THEN
*
*        Routines for general matrices
*
         CALL CTIMGE( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ),
     $                IWORK, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'GB' ) ) THEN
*
*        General band matrices
*
         IF( LDAMOK ) THEN
            CALL CTIMGB( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NNB, NBVAL,
     $                   NLDA, LDAVAL, TIMMIN, A( 1, 1 ),
     $                   A( LA-MKMAX+1, 3 ), IWORK, RESLTS, LDR1, LDR2,
     $                   LDR3, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )LINE( 1: 6 )
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'GT' ) ) THEN
*
*        Routines for general tridiagonal matrices
*
         CALL CTIMGT( LINE, NN, NVAL, NK, KVAL, NLDA, LDAVAL, TIMMIN,
     $                A( 1, 1 ), A( 1, 2 ), IWORK, RESLTS, LDR1, LDR2,
     $                LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'PO' ) ) THEN
*
*        Positive definite matrices
*
         CALL CTIMPO( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), IWORK,
     $                RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'PP' ) ) THEN
*
*        Positive definite packed matrices
*
         CALL CTIMPP( LINE, NN, NVAL, NK, KVAL, LA, TIMMIN, A( 1, 1 ),
     $                A( 1, 2 ), IWORK, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'PB' ) ) THEN
*
*        Positive definite banded matrices
*
         IF( LDAMOK ) THEN
            IF( MKMAX.GT.LA ) THEN
               J2 = 2
               I2 = 2*LA - MKMAX + 1
            ELSE
               J2 = 3
               I2 = LA - MKMAX + 1
            END IF
            CALL CTIMPB( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NNB, NBVAL,
     $                   NLDA, LDAVAL, TIMMIN, A( 1, 1 ), A( I2, J2 ),
     $                   IWORK, RESLTS, LDR1, LDR2, LDR3, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )LINE( 1: 6 )
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'PT' ) ) THEN
*
*        Routines for positive definite tridiagonal matrices
*
         CALL CTIMPT( LINE, NN, NVAL, NK, KVAL, NLDA, LDAVAL, TIMMIN, D,
     $                A( 1, 1 ), A( 1, 2 ), RESLTS, LDR1, LDR2, LDR3,
     $                NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'HE' ) ) THEN
*
*        Hermitian indefinite matrices
*
         CALL CTIMHE( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ),
     $                IWORK, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'HP' ) ) THEN
*
*        Hermitian indefinite packed matrices
*
         CALL CTIMHP( LINE, NN, NVAL, NK, KVAL, LA, TIMMIN, A( 1, 1 ),
     $                A( 1, 2 ), A( 1, 3 ), IWORK, RESLTS, LDR1, LDR2,
     $                LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'SY' ) ) THEN
*
*        Symmetric indefinite matrices
*
         CALL CTIMSY( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), A( 1, 3 ),
     $                IWORK, RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'SP' ) ) THEN
*
*        Symmetric indefinite packed matrices
*
         CALL CTIMSP( LINE, NN, NVAL, NK, KVAL, LA, TIMMIN, A( 1, 1 ),
     $                A( 1, 2 ), A( 1, 3 ), IWORK, RESLTS, LDR1, LDR2,
     $                LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'TR' ) ) THEN
*
*        Triangular matrices
*
         CALL CTIMTR( LINE, NN, NVAL, NK, KVAL, NNB, NBVAL, NLDA,
     $                LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ), RESLTS,
     $                LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'TP' ) ) THEN
*
*        Triangular packed matrices
*
         CALL CTIMTP( LINE, NN, NVAL, NK, KVAL, LA, TIMMIN, A( 1, 1 ),
     $                A( 1, 2 ), RESLTS, LDR1, LDR2, LDR3, NOUT )
*
      ELSE IF( LSAMEN( 2, C2, 'TB' ) ) THEN
*
*        Triangular band matrices
*
         IF( LDAMOK ) THEN
            IF( MKMAX.GT.LA ) THEN
               J2 = 2
               I2 = 2*LA - MKMAX + 1
            ELSE
               J2 = 3
               I2 = LA - MKMAX + 1
            END IF
            CALL CTIMTB( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NLDA,
     $                   LDAVAL, TIMMIN, A( 1, 1 ), A( I2, J2 ), RESLTS,
     $                   LDR1, LDR2, LDR3, NOUT )
         ELSE
            WRITE( NOUT, FMT = 9989 )LINE( 1: 6 )
         END IF
*
      ELSE IF( LSAMEN( 2, C2, 'LS' ) ) THEN
*
*        Least squares drivers
*
         CALL CTIMLS( LINE, NM, MVAL, NN, NVAL, NK, KVAL, NNB, NBVAL,
     $                NXVAL, NLDA, LDAVAL, TIMMIN, A( 1, 1 ), A( 1, 2 ),
     $                B( 1, 1 ), B( 1, 2 ), S, S( NMAX+1 ), OPCTBL,
     $                TIMTBL, FLPTBL, WORK, RWORK, IWORK, NOUT )
*
      ELSE
*
         WRITE( NOUT, FMT = 9987 )LINE( 1: 6 )
      END IF
*
*     Read the next line of the input file.
*
      READ( NIN, FMT = '(A)', END = 100 )LINE
      GO TO 90
*
*     Branch to this line when the last record is read.
*
  100 CONTINUE
      S2 = SECOND( )
      WRITE( NOUT, FMT = 9986 )
      WRITE( NOUT, FMT = 9985 )S2 - S1
  110 CONTINUE
*
 9999 FORMAT( ' Too many values of ', A, ' using ', A, ' = ', I2 )
 9998 FORMAT( ' *** LDA = ', I7, ' is too small, must have ',
     $      'LDA > 0.' )
 9997 FORMAT( ' *** ', A1, ' = ', I7, ' is too big:  ',
     $      'maximum allowed is', I7 )
 9996 FORMAT( ' *** N*NB is too big for N =', I6, ', NB =', I6,
     $      / ' --> Increase LA to at least ', I8 )
 9995 FORMAT( ' *** LDA*N is too big for the dense routines ', '(LDA =',
     $      I6, ', N =', I6, ')', / ' --> Increase LA to at least ',
     $      I8 )
 9994 FORMAT( ' *** (LDA+K)*M is too big for the band routines ',
     $      '(LDA=', I6, ', M=', I6, ', K=', I6, ')',
     $      / ' --> Increase LA to at least ', I8 )
 9993 FORMAT( ' The minimum time a subroutine will be timed = ', F6.3,
     $      ' seconds' )
 9992 FORMAT( ' The following parameter values will be used:' )
 9991 FORMAT( 4X, A7, 1X, 10I6, / 12X, 10I6 )
 9990 FORMAT( / ' ------------------------------',
     $      / ' >>>>>    Sample BLAS     <<<<<',
     $      / ' ------------------------------' )
 9989 FORMAT( 1X, A, ' not timed due to input errors', / )
 9988 FORMAT( / ' ------------------------------',
     $      / ' >>>>>    Timing data     <<<<<',
     $      / ' ------------------------------' )
 9987 FORMAT( 1X, A, ':  Unrecognized path or subroutine name', / )
 9986 FORMAT( ' End of tests' )
 9985 FORMAT( ' Total time used = ', F12.2, ' seconds' )
 9984 FORMAT( / ' Tests not done due to input errors' )
 9983 FORMAT( ' LAPACK VERSION 3.0, released June 30, 1999 ', / )
*
*     End of CTIMAA
*
      END
