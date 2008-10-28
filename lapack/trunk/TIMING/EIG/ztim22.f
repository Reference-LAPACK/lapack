      SUBROUTINE ZTIM22( LINE, NSIZES, NN, NTYPES, DOTYPE, NPARMS, NNB,
     $                   LDAS, TIMMIN, NOUT, ISEED, A, D, E, E2, U, URE,
     $                   UIM, TAU, TAURE, Z, ZRE, ZIM, WORK, LWORK,
     $                   RWORK, LLWORK, IWORK, TIMES, LDT1, LDT2, LDT3,
     $                   OPCNTS, LDO1, LDO2, LDO3, INFO )
*
*  -- LAPACK timing routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     October 2006
*
*     .. Scalar Arguments ..
      CHARACTER*80       LINE
      INTEGER            INFO, LDO1, LDO2, LDO3, LDT1, LDT2, LDT3,
     $                   LWORK, NOUT, NPARMS, NSIZES, NTYPES
      DOUBLE PRECISION   TIMMIN
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * ), LLWORK( * )
      INTEGER            ISEED( * ), IWORK( * ), LDAS( * ), NN( * ),
     $                   NNB( * )
      DOUBLE PRECISION   D( * ), E( * ), E2( * ),
     $                   OPCNTS( LDO1, LDO2, LDO3, * ), RWORK( * ),
     $                   TAURE( * ), TIMES( LDT1, LDT2, LDT3, * ),
     $                   UIM( * ), URE( * ), ZIM( * ), ZRE( * )
      COMPLEX*16         A( * ), TAU( * ), U( * ), WORK( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*     ZTIM22 times the LAPACK routines for the complex hermitian
*     eigenvalue problem.
*
*     For each N value in NN(1:NSIZES) and .TRUE. value in
*     DOTYPE(1:NTYPES), a matrix will be generated and used to test the
*     selected routines.  Thus, NSIZES*(number of .TRUE. values in
*     DOTYPE) matrices will be generated.
*
*  Arguments
*  =========
*
*  LINE    (input) CHARACTER*80
*          On entry, LINE contains the input line which requested
*          this routine.  This line may contain a subroutine name,
*          such as ZHETRD, indicating that only routine CHETRD will
*          be timed, or it may contain a generic name, such as ZST.
*          In this case, the rest of the line is scanned for the
*          first 12 non-blank characters, corresponding to the twelve
*          combinations of subroutine and options:
*          LAPACK:
*             1: ZHETRD
*             2: ZSTEQR(VECT='N')
*             3: ZUNGTR+ZSTEQR(VECT='V') (compare with IMTQL2+HTRIBK)
*             4: ZPTEQR(VECT='N')
*             5: ZUNGTR+ZPTEQR(VECT='V')
*             6. DSTEBZ+ZSTEIN+ZUNMTR
*             7. ZUNGTR+ZSTEDC(COMPQ='V')
*             8. ZSTEDC(COMPQ='I')+ZUNMTR
*             9. ZSTEGR(COMPQ='V')
*          EISPACK:
*            10: HTRIDI (compare with ZHETRD)
*            11: IMTQL1 (compare w/ ZSTEQR -- VECT='N')
*            12: IMTQL2+HTRIBK (compare w/ ZUNGTR+ZSTEQR(VECT='V') )
*          If a character is 'T' or 't', the corresponding routine in
*          this path is timed.  If the entire line is blank, all the
*          routines in the path are timed.
*
*  NSIZES  (input) INTEGER
*          The number of values of N contained in the vector NN.
*
*  NN      (input) INTEGER array, dimension( NSIZES )
*          The values of the matrix size N to be tested.  For each
*          N value in the array NN, and each .TRUE. value in DOTYPE,
*          a matrix A will be generated and used to test the routines.
*
*  NTYPES  (input) INTEGER
*          The number of types in DOTYPE.  Only the first MAXTYP
*          elements will be examined.  Exception: if NSIZES=1 and
*          NTYPES=MAXTYP+1, and DOTYPE=MAXTYP*f,t, then the input
*          value of A will be used.
*
*  DOTYPE  (input) LOGICAL
*          If DOTYPE(j) is .TRUE., then a matrix of type j will be
*          generated.  The matrix A has the form X**(-1) D X, where
*          X is unitary and D is diagonal with:
*          (j=1)  evenly spaced entries 1, ..., ULP with random signs.
*          (j=2)  geometrically spaced entries 1, ..., ULP with random
*                 signs.
*          (j=3)  "clustered" entries 1, ULP,..., ULP with random
*                 signs.
*          (j=4)  entries randomly chosen from ( ULP, 1 ).
*
*  NPARMS  (input) INTEGER
*          The number of values in each of the arrays NNB and LDAS.
*          For each matrix A generated according to NN and DOTYPE,
*          tests will be run with (NB,LDA)=
*          (NNB(1),LDAS(1)),...,(NNB(NPARMS), LDAS(NPARMS))
*
*  NNB     (input) INTEGER array, dimension( NPARMS )
*          The values of the blocksize ("NB") to be tested.
*
*  LDAS    (input) INTEGER array, dimension( NPARMS )
*          The values of LDA, the leading dimension of all matrices,
*          to be tested.
*
*  TIMMIN  (input) DOUBLE PRECISION
*          The minimum time a subroutine will be timed.
*
*  NOUT    (input) INTEGER
*          If NOUT > 0 then NOUT specifies the unit number
*          on which the output will be printed.  If NOUT <= 0, no
*          output is printed.
*
*  ISEED   (input/output) INTEGER array, dimension( 4 )
*          The random seed used by the random number generator, used
*          by the test matrix generator.  It is used and updated on
*          each call to ZTIM22
*
*  A       (workspace) COMPLEX*16 array, dimension( max(NN)*max(LDAS) )
*          The original matrix to be tested.
*
*  D       (workspace) DOUBLE PRECISION array, dimension( max(NN) )
*          The diagonal of the tridiagonal generated by ZHETRD/HTRIDI.
*
*  E       (workspace) DOUBLE PRECISION array, dimension( max(NN) )
*          The off-diagonal of the tridiagonal generated by
*          ZHETRD/HTRIDI.
*
*  E2      (workspace) DOUBLE PRECISION array, dimension( max(NN) )
*          The diagonal of a positive definite tridiagonal matrix
*          sent to ZPTEQR.  The off-diagonal is in array E.
*
*  U       (workspace) COMPLEX*16 array, dimension( max(NN)*max(LDAS) )
*          The array of Householder vectors output by ZHETRD.  This
*          array is used only when URE and UIM are not; thus, on
*          nearly all computers, URE may be EQUIVALENCEd with the
*          first half of U in the main (calling) routine, and UIM with
*          the second half, although this is a violation of the
*          FORTRAN-77 standard.
*
*  URE     (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          The array of the real parts of Householder vectors output by
*          HTRIDI.  This array is used only when U is not -- see the
*          note description of U.
*
*  UIM     (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          The array of the imaginary parts of Householder vectors
*          output by HTRIDI.  This array is used only when U is not --
*          see the description of U.
*
*  TAU     (workspace) COMPLEX*16 array, dimension( max(NN) )
*          The vector of coefficients for the Householder
*          transformations output by ZHETRD.  This array is used only
*          when TAURE is not; thus, on nearly all computers, TAURE may
*          be EQUIVALENCEd with TAU in the main (calling) routine,
*          although this is a violation of the FORTRAN-77 standard.
*
*  TAURE   (workspace) DOUBLE PRECISION array, dimension( 2*max(NN) )
*          The vector of complex (modulus 1) factors output by HTRIDI.
*          This vector is used only when TAU is not -- see the
*          description of TAU.
*
*  Z       (workspace) COMPLEX*16 array, dimension( max(NN)*max(LDAS) )
*          Various output arrays.  This array is used only when ZRE
*          and ZIM are not; thus, on nearly all computers, ZRE may be
*          EQUIVALENCEd with the first half of Z in the main (calling)
*          routine, and ZIM with the second half, although this is a
*          violation of the FORTRAN-77 standard.
*
*  ZRE     (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays (real parts).  This array is used
*          only when Z is not -- see the description of Z.
*
*  ZIM     (workspace) DOUBLE PRECISION array,
*                      dimension( max(NN)*max(LDAS) )
*          Various output arrays (imaginary parts).  This array is
*          used only when Z is not -- see the description of Z.
*
*  WORK    (workspace) COMPLEX*16 array, dimension( LWORK )
*
*  LWORK   (input) INTEGER
*          Number of elements in WORK.  It must be at least
*          max( (NNB + 2 )*LDAS, max(LDAS)*max(LDAS) )
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension
*                   ( max( 6*max(LDAS), NSIZES*NTYPES*NPARMS ),
*                     ( 1 + 3 * M + 2 * M * lg M + 3 * M**2 ) ),
*          where  M = max(lDAS), and lg M is the smallest integer k
*          such that 2^k >= N.
*          This should *not* be equivalenced to other arrays.
*
*  LLWORK  (workspace) LOGICAL array, dimension( NPARMS )
*
*  IWORK   (workspace) INTEGER array, dimension max( 5*max(LDAS),
*          ( 6 + 6*M + 5 * M * lg M ) ).
*
*  TIMES   (workspace) DOUBLE PRECISION array,
*                      dimension (LDT1,LDT2,LDT3,NSUBS)
*          TIMES(i,j,k,l) will be set to the run time (in seconds) for
*          subroutine l, with N=NN(k), matrix type j, and LDA=LDAS(i),
*          NBLOCK=NNB(i).
*
*  LDT1    (input) INTEGER
*          The first dimension of TIMES.  LDT1 >= min( 1, NPARMS ).
*
*  LDT2    (input) INTEGER
*          The second dimension of TIMES.  LDT2 >= min( 1, NTYPES ).
*
*  LDT3    (input) INTEGER
*          The third dimension of TIMES.  LDT3 >= min( 1, NSIZES ).
*
*  OPCNTS  (output) DOUBLE PRECISION array,
*                   dimension (LDO1,LDO2,LDO3,NSUBS)
*          OPCNTS(i,j,k,l) will be set to the number of floating-point
*          operations executed by subroutine l, with N=NN(k), matrix
*          type j, and LDA=LDAS(i), NBLOCK=NNB(i).
*
*  LDO1    (input) INTEGER
*          The first dimension of OPCNTS.  LDO1 >= min( 1, NPARMS ).
*
*  LDO2    (input) INTEGER
*          The second dimension of OPCNTS.  LDO2 >= min( 1, NTYPES ).
*
*  LDO3    (input) INTEGER
*          The third dimension of OPCNTS.  LDO3 >= min( 1, NSIZES ).
*
*  INFO    (output) INTEGER
*          Error flag.  It will be set to zero if no error occurred.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            MAXTYP, NSUBS
      PARAMETER          ( MAXTYP = 4, NSUBS = 12 )
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            RUNHTR, RUNTRD
      CHARACTER          UPLO
      INTEGER            I, IC, IINFO, IL, ILWORK, IMODE, IN, IPAR,
     $                   ISUB, ITYPE, IU, J, J1, J2, J3, J4, LASTL, LDA,
     $                   LDU, LGN, LIWEDC, LIWEVR, LRWEDC, LWEDC, LWEVR,
     $                   M, MTYPES, N, NB, NSPLIT, NANSOK, INFSOK
      DOUBLE PRECISION   ABSTOL, S1, S2, TIME, ULP, ULPINV, UNTIME, VL,
     $                   VU
*     ..
*     .. Local Arrays ..
      LOGICAL            TIMSUB( NSUBS )
      CHARACTER*4        PNAMES( 4 )
      CHARACTER*20       SUBNAM( NSUBS )
      INTEGER            IDUMMA( 1 ), INPARM( NSUBS ), IOLDSD( 4 ),
     $                   KMODE( MAXTYP )
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DOPLA, DSECND
      EXTERNAL           DLAMCH, DOPLA, DSECND, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           ATIMIN, DCOPY, DLASET, DPRTBE, DSTEBZ, HTRIBK,
     $                   HTRIDI, IMTQL1, IMTQL2, XLAENV, ZHETRD, ZLACPY,
     $                   ZLATMS, ZPTEQR, ZSTEDC, ZSTEGR, ZSTEIN, ZSTEQR,
     $                   ZUNGTR, ZUNMTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, INT, LOG, MAX, MIN
*     ..
*     .. Common blocks ..
      COMMON             / LATIME / OPS, ITCNT
*     ..
*     .. Scalars in Common ..
      DOUBLE PRECISION   ITCNT, OPS
*     ..
*     .. Data statements ..
      DATA               SUBNAM / 'ZHETRD', 'ZSTEQR(N)',
     $                   'ZUNGTR+ZSTEQR(V)', 'ZPTEQR(N)',
     $                   'ZUNGTR+ZPTEQR(V)', 'DSTEBZ+ZSTEIN+ZUNMTR',
     $                   'ZUNGTR+ZSTEDC(V)', 'ZSTEDC(I)+ZUNMTR',
     $                   'ZSTEGR(V)', 'HTRIDI', 'IMTQL1',
     $                   'IMTQL2+HTRIBK' /
      DATA               INPARM / 2, 1, 2, 1, 2, 2, 1, 1, 1, 1, 1, 1 /
      DATA               PNAMES / 'LDA', 'NB', 'bad1', 'bad2' /
      DATA               KMODE / 4, 3, 1, 5 /
*     ..
*     .. Executable Statements ..
*
*
*     Extract the timing request from the input line.
*
      CALL ATIMIN( 'ZST', LINE, NSUBS, SUBNAM, TIMSUB, NOUT, INFO )
*
*     Disable timing of ZSTEGR if we're non-IEEE-754 compliant.
*
      NANSOK = ILAENV( 10, 'ZSTEGR', ' ', 0, 0, 0, 0 )
      INFSOK = ILAENV( 11, 'ZSTEGR', ' ', 0, 0, 0, 0 )
      IF( NANSOK.NE.1 .OR. INFSOK.NE.1 )  THEN
         TIMSUB(9) = .FALSE.
      END IF
*
      IF( INFO.NE.0 )
     $   RETURN
*
*     Check that N <= LDA for the input values.
*
      DO 20 J2 = 1, NSIZES
         DO 10 J1 = 1, NPARMS
            IF( NN( J2 ).GT.LDAS( J1 ) ) THEN
               INFO = -8
               WRITE( NOUT, FMT = 9999 )LINE( 1: 6 )
 9999          FORMAT( 1X, A, ' timing run not attempted -- N > LDA',
     $               / )
               RETURN
            END IF
   10    CONTINUE
   20 CONTINUE
*
*     Check LWORK
*
      ILWORK = 0
      DO 30 J1 = 1, NPARMS
         ILWORK = MAX( ILWORK, ( NNB( J1 )+2 )*LDAS( J1 ) )
   30 CONTINUE
      IF( ILWORK.GT.LWORK ) THEN
         INFO = -18
         WRITE( NOUT, FMT = 9998 )LINE( 1: 6 )
 9998    FORMAT( 1X, A, ' timing run not attempted -- LWORK too small.',
     $         / )
         RETURN
      END IF
*
*     Check to see whether ZHETRD must be run.
*
*     RUNTRD -- if ZHETRD must be run.
*
      RUNTRD = .FALSE.
      IF( TIMSUB( 2 ) .OR. TIMSUB( 3 ) .OR. TIMSUB( 4 ) .OR.
     $    TIMSUB( 5 ) .OR. TIMSUB( 6 ) .OR. TIMSUB( 7 ) .OR.
     $    TIMSUB( 8 ) .OR. TIMSUB( 9 ) )RUNTRD = .TRUE.
*
*     Check to see whether HTRIDI must be run.
*
*     RUNHTR -- if HTRIDI must be run.
*
      RUNHTR = .FALSE.
      IF( TIMSUB( 10 ) .OR. TIMSUB( 11 ) .OR. TIMSUB( 12 ) )
     $   RUNHTR = .TRUE.
*
*     Various Constants
*
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      ULPINV = ONE / ULP
      CALL XLAENV( 9, 25 )
*
*     Zero out OPCNTS, TIMES
*
      DO 70 J4 = 1, NSUBS
         DO 60 J3 = 1, NSIZES
            DO 50 J2 = 1, NTYPES
               DO 40 J1 = 1, NPARMS
                  OPCNTS( J1, J2, J3, J4 ) = ZERO
                  TIMES( J1, J2, J3, J4 ) = ZERO
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
*
*     Do for each value of N:
*
      DO 650 IN = 1, NSIZES
*
         N = NN( IN )
         IF( N.GT.0 ) THEN
            LGN = INT( LOG( DBLE( N ) ) / LOG( TWO ) )
            IF( 2**LGN.LT.N )
     $         LGN = LGN + 1
            IF( 2**LGN.LT.N )
     $         LGN = LGN + 1
            LWEDC = 1 + 4*N + 2*N*LGN + 3*N**2
            LRWEDC = 1 + 3*N + 2*N*LGN + 3*N**2
            LIWEDC = 6 + 6*N + 5*N*LGN
            LWEVR = 18*N
            LIWEVR = 10*N
         ELSE
            LWEDC = 8
            LRWEDC = 7
            LIWEDC = 12
            LWEVR = 1
            LIWEVR = 1
         END IF
*
*        Do for each .TRUE. value in DOTYPE:
*
         MTYPES = MIN( MAXTYP, NTYPES )
         IF( NTYPES.EQ.MAXTYP+1 .AND. NSIZES.EQ.1 )
     $      MTYPES = NTYPES
         DO 640 ITYPE = 1, MTYPES
            IF( .NOT.DOTYPE( ITYPE ) )
     $         GO TO 640
*
*           Save random number seed for error messages
*
            DO 80 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   80       CONTINUE
*
*-----------------------------------------------------------------------
*
*           Time the LAPACK Routines
*
*           Generate A
*
            UPLO = 'L'
            IF( ITYPE.LE.MAXTYP ) THEN
               IMODE = KMODE( ITYPE )
               CALL ZLATMS( N, N, 'S', ISEED, 'S', RWORK, IMODE, ULPINV,
     $                      ONE, N, N, UPLO, A, N, WORK, IINFO )
               IF( IINFO.NE.0 ) THEN
                  WRITE( NOUT, FMT = 9997 )'ZLATMS', IINFO, N, ITYPE,
     $               0, IOLDSD
                  INFO = ABS( IINFO )
                  GO TO 640
               END IF
            END IF
*
*           Time ZHETRD for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 1 ) ) THEN
               DO 110 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time ZHETRD
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
   90             CONTINUE
                  CALL ZLACPY( UPLO, N, N, A, N, U, LDA )
                  CALL ZHETRD( UPLO, N, U, LDA, D, E, TAU, WORK, LWORK,
     $                         IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 190
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 90
*
*                 Subtract the time used in ZLACPY.
*
                  S1 = DSECND( )
                  DO 100 J = 1, IC
                     CALL ZLACPY( UPLO, N, N, A, N, Z, LDA )
  100             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 1 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 1 ) = DOPLA( 'ZHETRD', N, 0,
     $               0, 0, NB )
                  LDU = LDA
  110          CONTINUE
            ELSE
               IF( RUNTRD ) THEN
                  CALL ZLACPY( UPLO, N, N, A, N, U, N )
                  CALL ZHETRD( UPLO, N, U, N, D, E, TAU, WORK, LWORK,
     $                         IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 1 ), IINFO, N,
     $                  ITYPE, 0, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 190
                  END IF
                  LDU = N
               END IF
            END IF
*
*           Time ZSTEQR for each distinct LDA=LDAS(j)
*
            IF( TIMSUB( 2 ) ) THEN
               DO 150 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 120 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  120             CONTINUE
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time ZSTEQR with VECT='N'
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  130                CONTINUE
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL ZSTEQR( 'N', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                            RWORK( 2*LDA+1 ), IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 2 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 150
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 130
*
*                    Subtract the time used in DCOPY.
*
                     S1 = DSECND( )
                     DO 140 J = 1, IC
                        CALL DCOPY( N, D, 1, RWORK, 1 )
                        CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  140                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 2 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 2 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 2 )
                     TIMES( IPAR, ITYPE, IN, 2 ) = TIMES( LASTL, ITYPE,
     $                  IN, 2 )
                  END IF
  150          CONTINUE
            END IF
*
*           Time ZUNGTR + ZSTEQR(VECT='V') for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 3 ) ) THEN
               DO 180 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time ZUNGTR + ZSTEQR
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  160             CONTINUE
                  CALL ZLACPY( 'L', N, N, A, N, Z, LDA )
                  CALL ZUNGTR( 'L', N, Z, LDA, TAU, WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'ZUNGTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 180
                  END IF
                  CALL DCOPY( N, D, 1, RWORK, 1 )
                  CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL ZSTEQR( 'V', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                         RWORK( 2*LDA+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 3 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 180
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 160
*
*                 Subtract the time used in ZLACPY.
*
                  S1 = DSECND( )
                  DO 170 J = 1, IC
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL ZLACPY( 'L', N, N, A, N, Z, LDA )
  170             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 3 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 3 ) = OPS / DBLE( IC )
                  LDU = LDA
  180          CONTINUE
            END IF
*
  190       CONTINUE
*
*           Time ZPTEQR for each distinct LDA=LDAS(j)
*
            IF( TIMSUB( 4 ) ) THEN
               DO 240 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 200 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  200             CONTINUE
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time ZPTEQR with VECT='N'
*
*
*                    Modify the tridiagonal matrix to make it
*                    positive definite.
                     E2( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     DO 210 I = 2, N - 1
                        E2( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                            ABS( E( I-1 ) )
  210                CONTINUE
                     E2( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  220                CONTINUE
                     CALL DCOPY( N, E2, 1, RWORK( 1 ), 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL ZPTEQR( 'N', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                            RWORK( 2*LDA+1 ), IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 4 ), IINFO, N,
     $                     ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 240
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 220
*
*                    Subtract the time used in DCOPY.
*
                     S1 = DSECND( )
                     DO 230 J = 1, IC
                        CALL DCOPY( N, E2, 1, RWORK, 1 )
                        CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  230                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 4 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 4 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 4 )
                     TIMES( IPAR, ITYPE, IN, 4 ) = TIMES( LASTL, ITYPE,
     $                  IN, 4 )
                  END IF
  240          CONTINUE
            END IF
*
*           Time ZUNGTR + ZPTEQR(VECT='V') for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 5 ) ) THEN
               DO 290 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time ZUNGTR + ZPTEQR
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  250             CONTINUE
                  CALL ZLACPY( 'L', N, N, A, N, Z, LDA )
                  CALL ZUNGTR( 'L', N, Z, LDA, TAU, WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'ZUNGTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 290
                  END IF
*
*                 Modify the tridiagonal matrix to make it
*                 positive definite.
                  E2( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                  DO 260 I = 2, N - 1
                     E2( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                         ABS( E( I-1 ) )
  260             CONTINUE
                  E2( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
*
                  CALL DCOPY( N, E2, 1, RWORK, 1 )
                  CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL ZPTEQR( 'V', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                         RWORK( 2*LDA+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 5 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 290
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 250
*
*                 Subtract the time used in ZLACPY.
*
                  S1 = DSECND( )
                  DO 280 J = 1, IC
                     E2( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     DO 270 I = 2, N - 1
                        E2( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                            ABS( E( I-1 ) )
  270                CONTINUE
                     E2( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
*
                     CALL DCOPY( N, E2, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL ZLACPY( 'L', N, N, A, N, Z, LDA )
  280             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 5 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 5 ) = OPS / DBLE( IC )
                  LDU = LDA
  290          CONTINUE
            END IF
*
*           Time DSTEBZ+ZSTEIN+ZUNMTR for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 6 ) ) THEN
               VL = ZERO
               VU = ZERO
               IL = 1
               IU = N
               ABSTOL = ZERO
               DO 310 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
*
*                 Time DSTEBZ + ZSTEIN + ZUNMTR
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  300             CONTINUE
*
                  CALL DSTEBZ( 'A', 'B', N, VL, VU, IL, IU, ABSTOL, D,
     $                         E, M, NSPLIT, RWORK( 1 ), IWORK( 1 ),
     $                         IWORK( N+1 ), RWORK( 2*N+1 ),
     $                         IWORK( 2*N+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'DSTEBZ', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 310
                  END IF
*
                  CALL ZSTEIN( N, D, E, N, RWORK( 1 ), IWORK( 1 ),
     $                         IWORK( N+1 ), Z, LDA, RWORK( N+1 ),
     $                         IWORK( 2*N+1 ), IWORK( 3*N+1 ), IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 6 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 310
                  END IF
*
                  CALL ZUNMTR( 'L', 'L', 'N', N, N, U, LDU, TAU, Z, LDA,
     $                         WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'ZUNMTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 310
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 300
                  UNTIME = ZERO
*
                  TIMES( IPAR, ITYPE, IN, 6 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 6 ) = OPS / DBLE( IC )
                  LDU = LDA
  310          CONTINUE
            END IF
*
*           Time ZUNGTR + ZSTEDC(COMPQ='V') for each pair NNB(j),
*           LDAS(j)
*
            IF( TIMSUB( 7 ) ) THEN
               DO 340 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time ZUNGTR + ZSTEDC
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  320             CONTINUE
                  CALL ZLACPY( 'L', N, N, A, N, Z, LDA )
                  CALL ZUNGTR( 'L', N, Z, LDA, TAU, WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'ZUNGTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 340
                  END IF
                  CALL DCOPY( N, D, 1, RWORK, 1 )
                  CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL ZSTEDC( 'V', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                         WORK, LWEDC, RWORK( 2*LDA+1 ), LRWEDC,
     $                         IWORK, LIWEDC, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 7 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 340
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 320
*
*                 Subtract the time used in ZLACPY.
*
                  S1 = DSECND( )
                  DO 330 J = 1, IC
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL ZLACPY( 'L', N, N, A, N, Z, LDA )
  330             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 7 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 7 ) = OPS / DBLE( IC )
                  LDU = LDA
  340          CONTINUE
            END IF
*
*           Time ZSTEDC(COMPQ='I') + ZUNMTR for each pair NNB(j),
*           LDAS(j)
*
            IF( TIMSUB( 8 ) ) THEN
               DO 370 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
*                 Time  ZSTEDC + ZUNMTR
*
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  350             CONTINUE
                  CALL DCOPY( N, D, 1, RWORK, 1 )
                  CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL ZSTEDC( 'I', N, RWORK, RWORK( LDA+1 ), Z, LDA,
     $                         WORK, LWEDC, RWORK( 2*LDA+1 ), LRWEDC,
     $                         IWORK, LIWEDC, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 8 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 370
                  END IF
*
                  CALL ZUNMTR( 'L', 'L', 'N', N, N, U, LDU, TAU, Z, LDA,
     $                         WORK, LWORK, IINFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )'ZUNMTR', IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 370
                  END IF
*
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 350
*
*                 Subtract the time used in DCOPY.
*
                  S1 = DSECND( )
                  DO 360 J = 1, IC
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  360             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 8 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 8 ) = OPS / DBLE( IC )
                  LDU = LDA
  370          CONTINUE
            END IF
*
*           Time ZSTEGR(COMPQ='V') for each pair NNB(j), LDAS(j)
*
            IF( TIMSUB( 9 ) ) THEN
               DO 400 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
                  NB = MIN( N, NNB( IPAR ) )
                  CALL XLAENV( 1, NB )
                  CALL XLAENV( 2, 2 )
                  CALL XLAENV( 3, NB )
*
                  ABSTOL = ZERO
                  VL = ZERO
                  VU = ZERO
                  IL = 1
                  IU = N
                  IC = 0
                  OPS = ZERO
                  S1 = DSECND( )
  380             CONTINUE
                  CALL DCOPY( N, D, 1, RWORK, 1 )
                  CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                  CALL ZSTEGR( 'V', 'A', N, RWORK, RWORK( LDA+1 ), VL,
     $                         VU, IL, IU, ABSTOL, M, RWORK( 2*LDA+1 ),
     $                         Z, LDA, IWORK, RWORK( 3*LDA+1 ), LWEVR,
     $                         IWORK( 2*LDA+1 ), LIWEVR, INFO )
                  IF( IINFO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9997 )SUBNAM( 9 ), IINFO, N,
     $                  ITYPE, IPAR, IOLDSD
                     INFO = ABS( IINFO )
                     GO TO 400
                  END IF
                  S2 = DSECND( )
                  TIME = S2 - S1
                  IC = IC + 1
                  IF( TIME.LT.TIMMIN )
     $               GO TO 380
*
*                 Subtract the time used in DCOPY.
*
                  S1 = DSECND( )
                  DO 390 J = 1, IC
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  390             CONTINUE
                  S2 = DSECND( )
                  UNTIME = S2 - S1
*
                  TIMES( IPAR, ITYPE, IN, 9 ) = MAX( TIME-UNTIME,
     $               ZERO ) / DBLE( IC )
                  OPCNTS( IPAR, ITYPE, IN, 9 ) = OPS / DBLE( IC )
  400          CONTINUE
            END IF
*
*-----------------------------------------------------------------------
*
*           Time the EISPACK Routines
*
*           Skip routines if N <= 0 (EISPACK requirement)
*
            IF( N.LE.0 )
     $         GO TO 640
*
*           Time HTRIDI for each LDAS(j)
*
            IF( TIMSUB( 10 ) ) THEN
               DO 480 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 410 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  410             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time HTRIDI
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  420                CONTINUE
                     DO 440 J2 = 0, N - 1
                        DO 430 J1 = 1, N
                           URE( J1+LDA*J2 ) = DBLE( A( J1+N*J2 ) )
                           UIM( J1+LDA*J2 ) = DIMAG( A( J1+N*J2 ) )
  430                   CONTINUE
  440                CONTINUE
                     CALL HTRIDI( LDA, N, URE, UIM, D, E, RWORK, TAURE )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 420
*
*                    Subtract the time used in copying A.
*
                     S1 = DSECND( )
                     DO 470 J = 1, IC
                        DO 460 J2 = 0, N - 1
                           DO 450 J1 = 1, N
                              ZRE( J1+LDA*J2 ) = DBLE( A( J1+N*J2 ) )
                              ZIM( J1+LDA*J2 ) = DIMAG( A( J1+N*J2 ) )
  450                      CONTINUE
  460                   CONTINUE
  470                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
                     TIMES( IPAR, ITYPE, IN, 10 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 10 ) = OPS / DBLE( IC )
                     LDU = LDA
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 10 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 10 )
                     TIMES( IPAR, ITYPE, IN, 10 ) = TIMES( LASTL, ITYPE,
     $                  IN, 10 )
                  END IF
  480          CONTINUE
            ELSE
               IF( RUNHTR ) THEN
                  DO 500 J2 = 0, N - 1
                     DO 490 J1 = 1, N
                        URE( J1+N*J2 ) = DBLE( A( J1+N*J2 ) )
                        UIM( J1+N*J2 ) = DIMAG( A( J1+N*J2 ) )
  490                CONTINUE
  500             CONTINUE
                  CALL HTRIDI( N, N, URE, UIM, D, E, RWORK, TAURE )
                  LDU = N
               END IF
            END IF
*
*           Time IMTQL1 for each LDAS(j)
*
            IF( TIMSUB( 11 ) ) THEN
               DO 540 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 510 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  510             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Time IMTQL1
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  520                CONTINUE
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL IMTQL1( N, RWORK, RWORK( LDA+1 ), IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 11 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 550
                     END IF
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 520
*
*                    Subtract the time used in DCOPY
*
                     S1 = DSECND( )
                     DO 530 J = 1, IC
                        CALL DCOPY( N, D, 1, RWORK, 1 )
                        CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
  530                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 11 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 11 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 11 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 11 )
                     TIMES( IPAR, ITYPE, IN, 11 ) = TIMES( LASTL, ITYPE,
     $                  IN, 11 )
                  END IF
  540          CONTINUE
            END IF
  550       CONTINUE
*
*           Time IMTQL2 + HTRIBK for each LDAS(j)
*
            IF( TIMSUB( 12 ) ) THEN
               DO 630 IPAR = 1, NPARMS
                  LDA = LDAS( IPAR )
*
*                 If this value of LDA has come up before, just use
*                 the value previously computed.
*
                  LASTL = 0
                  DO 560 J = 1, IPAR - 1
                     IF( LDA.EQ.LDAS( J ) )
     $                  LASTL = J
  560             CONTINUE
*
                  IF( LASTL.EQ.0 ) THEN
*
*                    Change leading dimension of U
*
                     IF( LDA.GT.LDU ) THEN
                        DO 580 J2 = N - 1, 1, -1
                           DO 570 J1 = N, 1, -1
                              URE( J1+LDA*J2 ) = URE( J1+LDU*J2 )
                              UIM( J1+LDA*J2 ) = UIM( J1+LDU*J2 )
  570                      CONTINUE
  580                   CONTINUE
                        LDU = LDA
                     ELSE IF( LDA.LT.LDU ) THEN
                        DO 600 J2 = 1, N - 1
                           DO 590 J1 = 1, N
                              URE( J1+LDA*J2 ) = URE( J1+LDU*J2 )
                              UIM( J1+LDA*J2 ) = UIM( J1+LDU*J2 )
  590                      CONTINUE
  600                   CONTINUE
                        LDU = LDA
                     END IF
*
*                    Time IMTQL2 + HTRIBK
*
                     IC = 0
                     OPS = ZERO
                     S1 = DSECND( )
  610                CONTINUE
                     CALL DCOPY( N, D, 1, RWORK, 1 )
                     CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                     CALL DLASET( 'Full', N, N, ZERO, ONE, ZRE, LDA )
                     CALL IMTQL2( LDA, N, RWORK, RWORK( LDA+1 ), ZRE,
     $                            IINFO )
                     IF( IINFO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9997 )SUBNAM( 12 ), IINFO,
     $                     N, ITYPE, IPAR, IOLDSD
                        INFO = ABS( IINFO )
                        GO TO 640
                     END IF
                     CALL HTRIBK( LDA, N, URE, UIM, TAURE, N, ZRE, ZIM )
                     S2 = DSECND( )
                     TIME = S2 - S1
                     IC = IC + 1
                     IF( TIME.LT.TIMMIN )
     $                  GO TO 610
*
*                    Subtract the time used in copying
*
                     S1 = DSECND( )
                     DO 620 J = 1, IC
                        CALL DCOPY( N, D, 1, RWORK, 1 )
                        CALL DCOPY( N-1, E, 1, RWORK( LDA+1 ), 1 )
                        CALL DLASET( 'Full', N, N, ZERO, ONE, ZRE, LDA )
  620                CONTINUE
                     S2 = DSECND( )
                     UNTIME = S2 - S1
*
                     TIMES( IPAR, ITYPE, IN, 12 ) = MAX( TIME-UNTIME,
     $                  ZERO ) / DBLE( IC )
                     OPCNTS( IPAR, ITYPE, IN, 12 ) = OPS / DBLE( IC )
                  ELSE
                     OPCNTS( IPAR, ITYPE, IN, 12 ) = OPCNTS( LASTL,
     $                  ITYPE, IN, 12 )
                     TIMES( IPAR, ITYPE, IN, 12 ) = TIMES( LASTL, ITYPE,
     $                  IN, 12 )
                  END IF
  630          CONTINUE
            END IF
*
  640    CONTINUE
  650 CONTINUE
*
*-----------------------------------------------------------------------
*
*     Print a table of results for each timed routine.
*
      DO 660 ISUB = 1, NSUBS
         IF( TIMSUB( ISUB ) ) THEN
            CALL DPRTBE( SUBNAM( ISUB ), MTYPES, DOTYPE, NSIZES, NN,
     $                   INPARM( ISUB ), PNAMES, NPARMS, LDAS, NNB,
     $                   IDUMMA, IDUMMA, OPCNTS( 1, 1, 1, ISUB ), LDO1,
     $                   LDO2, TIMES( 1, 1, 1, ISUB ), LDT1, LDT2,
     $                   RWORK, LLWORK, NOUT )
         END IF
  660 CONTINUE
*
 9997 FORMAT( ' ZTIM22: ', A, ' returned INFO=', I6, '.', / 9X, 'N=',
     $      I6, ', ITYPE=', I6, ', IPAR=', I6, ', ISEED=(',
     $      3( I5, ',' ), I5, ')' )
*
      RETURN
*
*     End of ZTIM22
*
      END
