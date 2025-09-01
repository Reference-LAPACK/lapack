*> \brief \b ZCHKLQ
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL,
*                          NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AL, AC,
*                          B, X, XACT, TAU, WORK, RWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NM, NMAX, NN, NNB, NOUT, NRHS
*       DOUBLE PRECISION   THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            MVAL( * ), NBVAL( * ), NVAL( * ),
*      $                   NXVAL( * )
*       DOUBLE PRECISION   RWORK( * )
*       COMPLEX*16         A( * ), AC( * ), AF( * ), AL( * ), AQ( * ),
*      $                   B( * ), TAU( * ), WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZCHKLQ tests ZGELQF, ZUNGLQ and ZUNMLQ.
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
*> \param[in] NM
*> \verbatim
*>          NM is INTEGER
*>          The number of values of M contained in the vector MVAL.
*> \endverbatim
*>
*> \param[in] MVAL
*> \verbatim
*>          MVAL is INTEGER array, dimension (NM)
*>          The values of the matrix row dimension M.
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
*>          The values of the matrix column dimension N.
*> \endverbatim
*>
*> \param[in] NNB
*> \verbatim
*>          NNB is INTEGER
*>          The number of values of NB and NX contained in the
*>          vectors NBVAL and NXVAL.  The blocking parameters are used
*>          in pairs (NB,NX).
*> \endverbatim
*>
*> \param[in] NBVAL
*> \verbatim
*>          NBVAL is INTEGER array, dimension (NNB)
*>          The values of the blocksize NB.
*> \endverbatim
*>
*> \param[in] NXVAL
*> \verbatim
*>          NXVAL is INTEGER array, dimension (NNB)
*>          The values of the crossover point NX.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand side vectors to be generated for
*>          each linear system.
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
*>          The maximum value permitted for M or N, used in dimensioning
*>          the work arrays.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AQ
*> \verbatim
*>          AQ is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AL
*> \verbatim
*>          AL is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AC
*> \verbatim
*>          AC is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX*16 array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] XACT
*> \verbatim
*>          XACT is COMPLEX*16 array, dimension (NMAX*NRHS)
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16 array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX*16 array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension (NMAX)
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
      SUBROUTINE ZCHKLQ( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL,
     $                   NRHS, THRESH, TSTERR, NMAX, A, AF, AQ, AL, AC,
     $                   B, X, XACT, TAU, WORK, RWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NM, NMAX, NN, NNB, NOUT, NRHS
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            MVAL( * ), NBVAL( * ), NVAL( * ),
     $                   NXVAL( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( * ), AC( * ), AF( * ), AL( * ), AQ( * ),
     $                   B( * ), TAU( * ), WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 7 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 8 )
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          DIST, TYPE
      CHARACTER*3        PATH
      INTEGER            I, IK, IM, IMAT, IN, INB, INFO, K, KL, KU, LDA,
     $                   LWORK, M, MINMN, MODE, N, NB, NERRS, NFAIL, NK,
     $                   NRUN, NT, NX
      DOUBLE PRECISION   ANORM, CNDNUM
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 ), KVAL( 4 )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, XLAENV, ZERRLQ, ZGELS,
     $                   ZGET02, ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZLQT01,
     $                   ZLQT02, ZLQT03
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
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
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'LQ'
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
     $   CALL ZERRLQ( PATH, NOUT )
      INFOT = 0
      CALL XLAENV( 2, 2 )
*
      LDA = NMAX
      LWORK = NMAX*MAX( NMAX, NRHS )
*
*     Do for each value of M in MVAL.
*
      DO 70 IM = 1, NM
         M = MVAL( IM )
*
*        Do for each value of N in NVAL.
*
         DO 60 IN = 1, NN
            N = NVAL( IN )
            MINMN = MIN( M, N )
            DO 50 IMAT = 1, NTYPES
*
*              Do the tests only if DOTYPE( IMAT ) is true.
*
               IF( .NOT.DOTYPE( IMAT ) )
     $            GO TO 50
*
*              Set up parameters with ZLATB4 and generate a test matrix
*              with ZLATMS.
*
               CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
*
               SRNAMT = 'ZLATMS'
               CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, 'No packing', A, LDA,
     $                      WORK, INFO )
*
*              Check error code from ZLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 50
               END IF
*
*              Set some values for K: the first value must be MINMN,
*              corresponding to the call of ZLQT01; other values are
*              used in the calls of ZLQT02, and must not exceed MINMN.
*
               KVAL( 1 ) = MINMN
               KVAL( 2 ) = 0
               KVAL( 3 ) = 1
               KVAL( 4 ) = MINMN / 2
               IF( MINMN.EQ.0 ) THEN
                  NK = 1
               ELSE IF( MINMN.EQ.1 ) THEN
                  NK = 2
               ELSE IF( MINMN.LE.3 ) THEN
                  NK = 3
               ELSE
                  NK = 4
               END IF
*
*              Do for each value of K in KVAL
*
               DO 40 IK = 1, NK
                  K = KVAL( IK )
*
*                 Do for each pair of values (NB,NX) in NBVAL and NXVAL.
*
                  DO 30 INB = 1, NNB
                     NB = NBVAL( INB )
                     CALL XLAENV( 1, NB )
                     NX = NXVAL( INB )
                     CALL XLAENV( 3, NX )
                     DO I = 1, NTESTS
                        RESULT( I ) = ZERO
                     END DO
                     NT = 2
                     IF( IK.EQ.1 ) THEN
*
*                       Test ZGELQF
*
                        CALL ZLQT01( M, N, A, AF, AQ, AL, LDA, TAU,
     $                               WORK, LWORK, RWORK, RESULT( 1 ) )
                     ELSE IF( M.LE.N ) THEN
*
*                       Test ZUNGLQ, using factorization
*                       returned by ZLQT01
*
                        CALL ZLQT02( M, N, K, A, AF, AQ, AL, LDA, TAU,
     $                               WORK, LWORK, RWORK, RESULT( 1 ) )
                     END IF
                     IF( M.GE.K ) THEN
*
*                       Test ZUNMLQ, using factorization returned
*                       by ZLQT01
*
                        CALL ZLQT03( M, N, K, AF, AC, AL, AQ, LDA, TAU,
     $                               WORK, LWORK, RWORK, RESULT( 3 ) )
                        NT = NT + 4
*
*                       If M<=N and K=M, call ZGELS to solve a system
*                       with NRHS right hand sides and compute the
*                       residual.
*
                        IF( K.EQ.M .AND. INB.EQ.1 ) THEN
*
*                          Generate a solution and set the right
*                          hand side.
*
                           SRNAMT = 'ZLARHS'
                           CALL ZLARHS( PATH, 'New', 'Full',
     $                                  'No transpose', M, N, 0, 0,
     $                                  NRHS, A, LDA, XACT, LDA, B, LDA,
     $                                  ISEED, INFO )
*
                           CALL ZLACPY( 'Full', M, NRHS, B, LDA, X,
     $                                  LDA )
*
*                          Reset AF to the original matrix. ZGELS
*                          factors the matrix before solving the system.
*
                           CALL ZLACPY( 'Full', M, N, A, LDA, AF, LDA )
*
                           SRNAMT = 'ZGELS'
                           CALL ZGELS( 'No transpose', M, N, NRHS, AF,
     $                                 LDA, X, LDA, WORK, LWORK, INFO )
*
*                          Check error code from ZGELS.
*
                           IF( INFO.NE.0 )
     $                        CALL ALAERH( PATH, 'ZGELS', INFO, 0, 'N',
     $                                     M, N, NRHS, -1, NB, IMAT,
     $                                     NFAIL, NERRS, NOUT )
*
                           CALL ZGET02( 'No transpose', M, N, NRHS, A,
     $                                  LDA, X, LDA, B, LDA, RWORK,
     $                                  RESULT( 7 ) )
                           NT = NT + 1
                        END IF
                     END IF
*
*                    Print information about the tests that did not
*                    pass the threshold.
*
                     DO 20 I = 1, NT
                        IF( RESULT( I ).GE.THRESH ) THEN
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                        CALL ALAHD( NOUT, PATH )
                           WRITE( NOUT, FMT = 9999 )M, N, K, NB, NX,
     $                        IMAT, I, RESULT( I )
                           NFAIL = NFAIL + 1
                        END IF
   20                CONTINUE
                     NRUN = NRUN + NT
   30             CONTINUE
   40          CONTINUE
   50       CONTINUE
   60    CONTINUE
   70 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' M=', I5, ', N=', I5, ', K=', I5, ', NB=', I4, ', NX=',
     $      I5, ', type ', I2, ', test(', I2, ')=', G12.5 )
      RETURN
*
*     End of ZCHKLQ
*
      END
