*> \brief \b ZCHKHE_RK
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZCHKHE_RK( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
*                             THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
*                             XACT, WORK, RWORK, IWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NMAX, NN, NNB, NNS, NOUT
*       DOUBLE PRECISION   THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), E( * ),
*      $                   WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZCHKHE_RK tests ZHETRF_RK, -TRI_3, -TRS_3,
*> and -CON_3.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DOTYPE
*> \verbatim
*>          DOTYPE is LOGICAL array, dimension (NTYPES)
*>          The matrix types to be used for testing.  Matrices of type j
*>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
*>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER
*>          The number of values of N contained in the vector NVAL.
*> \endverbatim
*>
*> \param[in] NVAL
*> \verbatim
*>          NVAL is INTEGER array, dimension (NN)
*>          The values of the matrix dimension N.
*> \endverbatim
*>
*> \param[in] NNB
*> \verbatim
*>          NNB is INTEGER
*>          The number of values of NB contained in the vector NBVAL.
*> \endverbatim
*>
*> \param[in] NBVAL
*> \verbatim
*>          NBVAL is INTEGER array, dimension (NNB)
*>          The values of the blocksize NB.
*> \endverbatim
*>
*> \param[in] NNS
*> \verbatim
*>          NNS is INTEGER
*>          The number of values of NRHS contained in the vector NSVAL.
*> \endverbatim
*>
*> \param[in] NSVAL
*> \verbatim
*>          NSVAL is INTEGER array, dimension (NNS)
*>          The values of the number of right hand sides NRHS.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is DOUBLE PRECISION
*>          The threshold value for the test ratios.  A result is
*>          included in the output file if RESULT >= THRESH.  To have
*>          every test ratio printed, use THRESH = 0.
*> \endverbatim
*>
*> \param[in] TSTERR
*> \verbatim
*>          TSTERR is LOGICAL
*>          Flag that indicates whether error exits are to be tested.
*> \endverbatim
*>
*> \param[in] NMAX
*> \verbatim
*>          NMAX is INTEGER
*>          The maximum value permitted for N, used in dimensioning the
*>          work arrays.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is CCOMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AFAC
*> \verbatim
*>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is COMPLEX*16 array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] AINV
*> \verbatim
*>          AINV is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is CCOMPLEX*16 array, dimension (NMAX*NSMAX)
*>          where NSMAX is the largest entry in NSVAL.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] XACT
*> \verbatim
*>          XACT is COMPLEX*16 array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (NMAX*max(3,NSMAX))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (max(NMAX,2*NSMAX)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (2*NMAX)
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The unit number for output.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup complex16_lin
*
*  =====================================================================
      SUBROUTINE ZCHKHE_RK( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
     $                      THRESH, TSTERR, NMAX, A, AFAC, E, AINV, B,
     $                      X, XACT, WORK, RWORK, IWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NMAX, NN, NNB, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), E( * ),
     $                   WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   ONEHALF
      PARAMETER          ( ONEHALF = 0.5D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 10 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 7 )
*     ..
*     .. Local Scalars ..
      LOGICAL            TRFCON, ZEROT
      CHARACTER          DIST, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH, MATPATH
      INTEGER            I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS,
     $                   ITEMP, ITEMP2, IUPLO, IZERO, J, K, KL, KU, LDA,
     $                   LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS,
     $                   NRUN, NT
      DOUBLE PRECISION   ALPHA, ANORM, CNDNUM, CONST, SING_MAX,
     $                   SING_MIN, RCOND, RCONDC, DTEMP
*     ..
*     .. Local Arrays ..
      CHARACTER          UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 ), IDUMMY( 1 )
      DOUBLE PRECISION   RESULT( NTESTS )
      COMPLEX*16         BLOCK( 2, 2 ), ZDUMMY( 1 )
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DGET06, ZLANGE, ZLANHE
      EXTERNAL           DGET06, ZLANGE, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, ZERRHE, ZGESVD, ZGET04,
     $                   ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZPOT02, ZPOT03,
     $                   ZHECON_3, ZHET01_3, ZHETRF_RK, ZHETRI_3,
     $                   ZHETRS_3, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN, SQRT
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      CHARACTER*32       SRNAMT
      INTEGER            INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
*
*     Test path
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'HK'
*
*     Path to generate matrices
*
      MATPATH( 1: 1 ) = 'Zomplex precision'
      MATPATH( 2: 3 ) = 'HE'
*
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
*     Test the error exits
*
      IF( TSTERR )
     $   CALL ZERRHE( PATH, NOUT )
      INFOT = 0
*
*     Set the minimum block size for which the block routine should
*     be used, which will be later returned by ILAENV
*
      CALL XLAENV( 2, 2 )
*
*     Do for each value of N in NVAL
*
      DO 270 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         IZERO = 0
*
*        Do for each value of matrix type IMAT
*
         DO 260 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 260
*
*           Skip types 3, 4, 5, or 6 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 260
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 250 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
*                 Begin generate the test matrix A.
*
*                 Set up parameters with ZLATB4 for the matrix generator
*                 based on the type of matrix to be generated.
*
                  CALL ZLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM,
     $                         MODE, CNDNUM, DIST )
*
*                 Generate a matrix with ZLATMS.
*
                  SRNAMT = 'ZLATMS'
                  CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                         CNDNUM, ANORM, KL, KU, UPLO, A, LDA,
     $                         WORK, INFO )
*
*                 Check error code from ZLATMS and handle error.
*
                  IF( INFO.NE.0 ) THEN
                     CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N,
     $                            -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
*
*                    Skip all tests for this generated matrix
*
                     GO TO 250
                  END IF
*
*                 For matrix types 3-6, zero one or more rows and
*                 columns of the matrix to test that INFO is returned
*                 correctly.
*
                  IF( ZEROT ) THEN
                     IF( IMAT.EQ.3 ) THEN
                        IZERO = 1
                     ELSE IF( IMAT.EQ.4 ) THEN
                        IZERO = N
                     ELSE
                        IZERO = N / 2 + 1
                     END IF
*
                     IF( IMAT.LT.6 ) THEN
*
*                       Set row and column IZERO to zero.
*
                        IF( IUPLO.EQ.1 ) THEN
                           IOFF = ( IZERO-1 )*LDA
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = CZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = CZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        ELSE
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = CZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = CZERO
   50                      CONTINUE
                        END IF
                     ELSE
                        IF( IUPLO.EQ.1 ) THEN
*
*                          Set the first IZERO rows and columns to zero.
*
                           IOFF = 0
                           DO 70 J = 1, N
                              I2 = MIN( J, IZERO )
                              DO 60 I = 1, I2
                                 A( IOFF+I ) = CZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                        ELSE
*
*                          Set the last IZERO rows and columns to zero.
*
                           IOFF = 0
                           DO 90 J = 1, N
                              I1 = MAX( J, IZERO )
                              DO 80 I = I1, N
                                 A( IOFF+I ) = CZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        END IF
                     END IF
                  ELSE
                     IZERO = 0
                  END IF
*
*                 End generate the test matrix A.
*
*
*              Do for each value of NB in NBVAL
*
               DO 240 INB = 1, NNB
*
*                 Set the optimal blocksize, which will be later
*                 returned by ILAENV.
*
                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
*
*                 Copy the test matrix A into matrix AFAC which
*                 will be factorized in place. This is needed to
*                 preserve the test matrix A for subsequent tests.
*
                  CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
*
*                 Compute the L*D*L**T or U*D*U**T factorization of the
*                 matrix. IWORK stores details of the interchanges and
*                 the block structure of D. AINV is a work array for
*                 block factorization, LWORK is the length of AINV.
*
                  LWORK = MAX( 2, NB )*LDA
                  SRNAMT = 'ZHETRF_RK'
                  CALL ZHETRF_RK( UPLO, N, AFAC, LDA, E, IWORK, AINV,
     $                            LWORK, INFO )
*
*                 Adjust the expected value of INFO to account for
*                 pivoting.
*
                  K = IZERO
                  IF( K.GT.0 ) THEN
  100                CONTINUE
                     IF( IWORK( K ).LT.0 ) THEN
                        IF( IWORK( K ).NE.-K ) THEN
                           K = -IWORK( K )
                           GO TO 100
                        END IF
                     ELSE IF( IWORK( K ).NE.K ) THEN
                        K = IWORK( K )
                        GO TO 100
                     END IF
                  END IF
*
*                 Check error code from ZHETRF_RK and handle error.
*
                  IF( INFO.NE.K)
     $               CALL ALAERH( PATH, 'ZHETRF_RK', INFO, K,
     $                            UPLO, N, N, -1, -1, NB, IMAT,
     $                            NFAIL, NERRS, NOUT )
*
*                 Set the condition estimate flag if the INFO is not 0.
*
                  IF( INFO.NE.0 ) THEN
                     TRFCON = .TRUE.
                  ELSE
                     TRFCON = .FALSE.
                  END IF
*
*+    TEST 1
*                 Reconstruct matrix from factors and compute residual.
*
                  CALL ZHET01_3( UPLO, N, A, LDA, AFAC, LDA, E, IWORK,
     $                           AINV, LDA, RWORK, RESULT( 1 ) )
                  NT = 1
*
*+    TEST 2
*                 Form the inverse and compute the residual,
*                 if the factorization was competed without INFO > 0
*                 (i.e. there is no zero rows and columns).
*                 Do it only for the first block size.
*
                  IF( INB.EQ.1 .AND. .NOT.TRFCON ) THEN
                     CALL ZLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
                     SRNAMT = 'ZHETRI_3'
*
*                    Another reason that we need to compute the inverse
*                    is that ZPOT03 produces RCONDC which is used later
*                    in TEST6 and TEST7.
*
                     LWORK = (N+NB+1)*(NB+3)
                     CALL ZHETRI_3( UPLO, N, AINV, LDA, E, IWORK, WORK,
     $                              LWORK, INFO )
*
*                    Check error code from ZHETRI_3 and handle error.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'ZHETRI_3', INFO, -1,
     $                               UPLO, N, N, -1, -1, -1, IMAT,
     $                               NFAIL, NERRS, NOUT )
*
*                    Compute the residual for a Hermitian matrix times
*                    its inverse.
*
                     CALL ZPOT03( UPLO, N, A, LDA, AINV, LDA, WORK, LDA,
     $                            RWORK, RCONDC, RESULT( 2 ) )
                     NT = 2
                  END IF
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 110 K = 1, NT
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K,
     $                     RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  110             CONTINUE
                  NRUN = NRUN + NT
*
*+    TEST 3
*                 Compute largest element in U or L
*
                  RESULT( 3 ) = ZERO
                  DTEMP = ZERO
*
                  CONST = ( ( ALPHA**2-ONE ) / ( ALPHA**2-ONEHALF ) ) /
     $                    ( ONE-ALPHA )
*
                  IF( IUPLO.EQ.1 ) THEN
*
*                 Compute largest element in U
*
                     K = N
  120                CONTINUE
                     IF( K.LE.1 )
     $                  GO TO 130
*
                     IF( IWORK( K ).GT.ZERO ) THEN
*
*                       Get max absolute value from elements
*                       in column k in U
*
                        DTEMP = ZLANGE( 'M', K-1, 1,
     $                          AFAC( ( K-1 )*LDA+1 ), LDA, RWORK )
                     ELSE
*
*                       Get max absolute value from elements
*                       in columns k and k-1 in U
*
                        DTEMP = ZLANGE( 'M', K-2, 2,
     $                          AFAC( ( K-2 )*LDA+1 ), LDA, RWORK )
                        K = K - 1
*
                     END IF
*
*                    DTEMP should be bounded by CONST
*
                     DTEMP = DTEMP - CONST + THRESH
                     IF( DTEMP.GT.RESULT( 3 ) )
     $                  RESULT( 3 ) = DTEMP
*
                     K = K - 1
*
                     GO TO 120
  130                CONTINUE
*
                  ELSE
*
*                 Compute largest element in L
*
                     K = 1
  140                CONTINUE
                     IF( K.GE.N )
     $                  GO TO 150
*
                     IF( IWORK( K ).GT.ZERO ) THEN
*
*                       Get max absolute value from elements
*                       in column k in L
*
                        DTEMP = ZLANGE( 'M', N-K, 1,
     $                          AFAC( ( K-1 )*LDA+K+1 ), LDA, RWORK )
                     ELSE
*
*                       Get max absolute value from elements
*                       in columns k and k+1 in L
*
                        DTEMP = ZLANGE( 'M', N-K-1, 2,
     $                          AFAC( ( K-1 )*LDA+K+2 ), LDA, RWORK )
                        K = K + 1
*
                     END IF
*
*                    DTEMP should be bounded by CONST
*
                     DTEMP = DTEMP - CONST + THRESH
                     IF( DTEMP.GT.RESULT( 3 ) )
     $                  RESULT( 3 ) = DTEMP
*
                     K = K + 1
*
                     GO TO 140
  150                CONTINUE
                  END IF
*
*
*+    TEST 4
*                 Compute largest 2-Norm (condition number)
*                 of 2-by-2 diag blocks
*
                  RESULT( 4 ) = ZERO
                  DTEMP = ZERO
*
                  CONST = ( ( ALPHA**2-ONE ) / ( ALPHA**2-ONEHALF ) )*
     $                    ( ( ONE + ALPHA ) / ( ONE - ALPHA ) )
                  CALL ZLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
*
                  IF( IUPLO.EQ.1 ) THEN
*
*                    Loop backward for UPLO = 'U'
*
                     K = N
  160                CONTINUE
                     IF( K.LE.1 )
     $                  GO TO 170
*
                     IF( IWORK( K ).LT.ZERO ) THEN
*
*                       Get the two singular values
*                       (real and non-negative) of a 2-by-2 block,
*                       store them in RWORK array
*
                        BLOCK( 1, 1 ) = AFAC( ( K-2 )*LDA+K-1 )
                        BLOCK( 1, 2 ) = E( K )
                        BLOCK( 2, 1 ) = DCONJG( BLOCK( 1, 2 ) )
                        BLOCK( 2, 2 ) = AFAC( (K-1)*LDA+K )
*
                        CALL ZGESVD( 'N', 'N', 2, 2, BLOCK, 2, RWORK,
     $                               ZDUMMY, 1, ZDUMMY, 1,
     $                               WORK, 6, RWORK( 3 ), INFO )
*
*
                        SING_MAX = RWORK( 1 )
                        SING_MIN = RWORK( 2 )
*
                        DTEMP = SING_MAX / SING_MIN
*
*                       DTEMP should be bounded by CONST
*
                        DTEMP = DTEMP - CONST + THRESH
                        IF( DTEMP.GT.RESULT( 4 ) )
     $                     RESULT( 4 ) = DTEMP
                        K = K - 1
*
                     END IF
*
                     K = K - 1
*
                     GO TO 160
  170                CONTINUE
*
                  ELSE
*
*                    Loop forward for UPLO = 'L'
*
                     K = 1
  180                CONTINUE
                     IF( K.GE.N )
     $                  GO TO 190
*
                     IF( IWORK( K ).LT.ZERO ) THEN
*
*                       Get the two singular values
*                       (real and non-negative) of a 2-by-2 block,
*                       store them in RWORK array
*
                        BLOCK( 1, 1 ) = AFAC( ( K-1 )*LDA+K )
                        BLOCK( 2, 1 ) = E( K )
                        BLOCK( 1, 2 ) = DCONJG( BLOCK( 2, 1 ) )
                        BLOCK( 2, 2 ) = AFAC( K*LDA+K+1 )
*
                        CALL ZGESVD( 'N', 'N', 2, 2, BLOCK, 2, RWORK,
     $                               ZDUMMY, 1, ZDUMMY, 1,
     $                               WORK, 6, RWORK(3), INFO )
*
                        SING_MAX = RWORK( 1 )
                        SING_MIN = RWORK( 2 )
*
                        DTEMP = SING_MAX / SING_MIN
*
*                       DTEMP should be bounded by CONST
*
                        DTEMP = DTEMP - CONST + THRESH
                        IF( DTEMP.GT.RESULT( 4 ) )
     $                     RESULT( 4 ) = DTEMP
                        K = K + 1
*
                     END IF
*
                     K = K + 1
*
                     GO TO 180
  190                CONTINUE
                  END IF
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  DO 200 K = 3, 4
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                     CALL ALAHD( NOUT, PATH )
                        WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K,
     $                     RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
  200             CONTINUE
                  NRUN = NRUN + 2
*
*                 Skip the other tests if this is not the first block
*                 size.
*
                  IF( INB.GT.1 )
     $               GO TO 240
*
*                 Do only the condition estimate if INFO is not 0.
*
                  IF( TRFCON ) THEN
                     RCONDC = ZERO
                     GO TO 230
                  END IF
*
*                 Do for each value of NRHS in NSVAL.
*
                  DO 220 IRHS = 1, NNS
                     NRHS = NSVAL( IRHS )
*
*                    Begin loop over NRHS values
*
*
*+    TEST 5 ( Using TRS_3)
*                 Solve and compute residual for  A * X = B.
*
*                    Choose a set of NRHS random solution vectors
*                    stored in XACT and set up the right hand side B
*
                     SRNAMT = 'ZLARHS'
                     CALL ZLARHS( MATPATH, XTYPE, UPLO, ' ', N, N,
     $                            KL, KU, NRHS, A, LDA, XACT, LDA,
     $                            B, LDA, ISEED, INFO )
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
*
                     SRNAMT = 'ZHETRS_3'
                     CALL ZHETRS_3( UPLO, N, NRHS, AFAC, LDA, E, IWORK,
     $                              X, LDA, INFO )
*
*                    Check error code from ZHETRS_3 and handle error.
*
                     IF( INFO.NE.0 )
     $                  CALL ALAERH( PATH, 'ZHETRS_3', INFO, 0,
     $                               UPLO, N, N, -1, -1, NRHS, IMAT,
     $                               NFAIL, NERRS, NOUT )
*
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
*
*                    Compute the residual for the solution
*
                     CALL ZPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $                            LDA, RWORK, RESULT( 5 ) )
*
*+    TEST 6
*                 Check solution from generated exact solution.
*
                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                            RESULT( 6 ) )
*
*                    Print information about the tests that did not pass
*                    the threshold.
*
                     DO 210 K = 5, 6
                        IF( RESULT( K ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS,
     $                        IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        END IF
  210                CONTINUE
                     NRUN = NRUN + 2
*
*                 End do for each value of NRHS in NSVAL.
*
  220             CONTINUE
*
*+    TEST 7
*                 Get an estimate of RCOND = 1/CNDNUM.
*
  230             CONTINUE
                  ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )
                  SRNAMT = 'ZHECON_3'
                  CALL ZHECON_3( UPLO, N, AFAC, LDA, E, IWORK, ANORM,
     $                           RCOND, WORK, INFO )
*
*                 Check error code from ZHECON_3 and handle error.
*
                  IF( INFO.NE.0 )
     $               CALL ALAERH( PATH, 'ZHECON_3', INFO, 0,
     $                            UPLO, N, N, -1, -1, -1, IMAT,
     $                            NFAIL, NERRS, NOUT )
*
*                 Compute the test ratio to compare values of RCOND
*
                  RESULT( 7 ) = DGET06( RCOND, RCONDC )
*
*                 Print information about the tests that did not pass
*                 the threshold.
*
                  IF( RESULT( 7 ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                  CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9997 )UPLO, N, IMAT, 7,
     $                  RESULT( 7 )
                     NFAIL = NFAIL + 1
                  END IF
                  NRUN = NRUN + 1
  240          CONTINUE
*
  250       CONTINUE
  260    CONTINUE
  270 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NB =', I4, ', type ',
     $      I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test ', I2, ', ratio =', G12.5 )
 9997 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ',', 10X, ' type ', I2,
     $      ', test ', I2, ', ratio =', G12.5 )
      RETURN
*
*     End of ZCHKHE_RK
*
      END
