*> \brief \b SCHKPT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
*                          A, D, E, B, X, XACT, WORK, RWORK, NOUT )
*
*       .. Scalar Arguments ..
*       LOGICAL            TSTERR
*       INTEGER            NN, NNS, NOUT
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            NSVAL( * ), NVAL( * )
*       REAL               A( * ), B( * ), D( * ), E( * ), RWORK( * ),
*      $                   WORK( * ), X( * ), XACT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SCHKPT tests SPTTRF, -TRS, -RFS, and -CON
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
*>          THRESH is REAL
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
*> \param[out] A
*> \verbatim
*>          A is REAL array, dimension (NMAX*2)
*> \endverbatim
*>
*> \param[out] D
*> \verbatim
*>          D is REAL array, dimension (NMAX*2)
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is REAL array, dimension (NMAX*2)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is REAL array, dimension (NMAX*NSMAX)
*>          where NSMAX is the largest entry in NSVAL.
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is REAL array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] XACT
*> \verbatim
*>          XACT is REAL array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension
*>                      (NMAX*max(3,NSMAX))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension
*>                      (max(NMAX,2*NSMAX))
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
*> \ingroup single_lin
*
*  =====================================================================
      SUBROUTINE SCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
     $                   A, D, E, B, X, XACT, WORK, RWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      INTEGER            NN, NNS, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            NSVAL( * ), NVAL( * )
      REAL               A( * ), B( * ), D( * ), E( * ), RWORK( * ),
     $                   WORK( * ), X( * ), XACT( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 12 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 7 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, TYPE
      CHARACTER*3        PATH
      INTEGER            I, IA, IMAT, IN, INFO, IRHS, IX, IZERO, J, K,
     $                   KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT,
     $                   NRHS, NRUN
      REAL               AINVNM, ANORM, COND, DMAX, RCOND, RCONDC
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      REAL               RESULT( NTESTS ), Z( 3 )
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SASUM, SGET06, SLANST
      EXTERNAL           ISAMAX, SASUM, SGET06, SLANST
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, SCOPY, SERRGT, SGET04,
     $                   SLACPY, SLAPTM, SLARNV, SLATB4, SLATMS, SPTCON,
     $                   SPTRFS, SPTT01, SPTT02, SPTT05, SPTTRF, SPTTRS,
     $                   SSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
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
      DATA               ISEEDY / 0, 0, 0, 1 /
*     ..
*     .. Executable Statements ..
*
      PATH( 1: 1 ) = 'Single precision'
      PATH( 2: 3 ) = 'PT'
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
     $   CALL SERRGT( PATH, NOUT )
      INFOT = 0
*
      DO 110 IN = 1, NN
*
*        Do for each value of N in NVAL.
*
         N = NVAL( IN )
         LDA = MAX( 1, N )
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 100 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( N.GT.0 .AND. .NOT.DOTYPE( IMAT ) )
     $         GO TO 100
*
*           Set up parameters with SLATB4.
*
            CALL SLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                   COND, DIST )
*
            ZEROT = IMAT.GE.8 .AND. IMAT.LE.10
            IF( IMAT.LE.6 ) THEN
*
*              Type 1-6:  generate a symmetric tridiagonal matrix of
*              known condition number in lower triangular band storage.
*
               SRNAMT = 'SLATMS'
               CALL SLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, COND,
     $                      ANORM, KL, KU, 'B', A, 2, WORK, INFO )
*
*              Check the error code from SLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'SLATMS', INFO, 0, ' ', N, N, KL,
     $                         KU, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               END IF
               IZERO = 0
*
*              Copy the matrix to D and E.
*
               IA = 1
               DO 20 I = 1, N - 1
                  D( I ) = A( IA )
                  E( I ) = A( IA+1 )
                  IA = IA + 2
   20          CONTINUE
               IF( N.GT.0 )
     $            D( N ) = A( IA )
            ELSE
*
*              Type 7-12:  generate a diagonally dominant matrix with
*              unknown condition number in the vectors D and E.
*
               IF( .NOT.ZEROT .OR. .NOT.DOTYPE( 7 ) ) THEN
*
*                 Let D and E have values from [-1,1].
*
                  CALL SLARNV( 2, ISEED, N, D )
                  CALL SLARNV( 2, ISEED, N-1, E )
*
*                 Make the tridiagonal matrix diagonally dominant.
*
                  IF( N.EQ.1 ) THEN
                     D( 1 ) = ABS( D( 1 ) )
                  ELSE
                     D( 1 ) = ABS( D( 1 ) ) + ABS( E( 1 ) )
                     D( N ) = ABS( D( N ) ) + ABS( E( N-1 ) )
                     DO 30 I = 2, N - 1
                        D( I ) = ABS( D( I ) ) + ABS( E( I ) ) +
     $                           ABS( E( I-1 ) )
   30                CONTINUE
                  END IF
*
*                 Scale D and E so the maximum element is ANORM.
*
                  IX = ISAMAX( N, D, 1 )
                  DMAX = D( IX )
                  CALL SSCAL( N, ANORM / DMAX, D, 1 )
                  CALL SSCAL( N-1, ANORM / DMAX, E, 1 )
*
               ELSE IF( IZERO.GT.0 ) THEN
*
*                 Reuse the last matrix by copying back the zeroed out
*                 elements.
*
                  IF( IZERO.EQ.1 ) THEN
                     D( 1 ) = Z( 2 )
                     IF( N.GT.1 )
     $                  E( 1 ) = Z( 3 )
                  ELSE IF( IZERO.EQ.N ) THEN
                     E( N-1 ) = Z( 1 )
                     D( N ) = Z( 2 )
                  ELSE
                     E( IZERO-1 ) = Z( 1 )
                     D( IZERO ) = Z( 2 )
                     E( IZERO ) = Z( 3 )
                  END IF
               END IF
*
*              For types 8-10, set one row and column of the matrix to
*              zero.
*
               IZERO = 0
               IF( IMAT.EQ.8 ) THEN
                  IZERO = 1
                  Z( 2 ) = D( 1 )
                  D( 1 ) = ZERO
                  IF( N.GT.1 ) THEN
                     Z( 3 ) = E( 1 )
                     E( 1 ) = ZERO
                  END IF
               ELSE IF( IMAT.EQ.9 ) THEN
                  IZERO = N
                  IF( N.GT.1 ) THEN
                     Z( 1 ) = E( N-1 )
                     E( N-1 ) = ZERO
                  END IF
                  Z( 2 ) = D( N )
                  D( N ) = ZERO
               ELSE IF( IMAT.EQ.10 ) THEN
                  IZERO = ( N+1 ) / 2
                  IF( IZERO.GT.1 ) THEN
                     Z( 1 ) = E( IZERO-1 )
                     E( IZERO-1 ) = ZERO
                     Z( 3 ) = E( IZERO )
                     E( IZERO ) = ZERO
                  END IF
                  Z( 2 ) = D( IZERO )
                  D( IZERO ) = ZERO
               END IF
            END IF
*
            CALL SCOPY( N, D, 1, D( N+1 ), 1 )
            IF( N.GT.1 )
     $         CALL SCOPY( N-1, E, 1, E( N+1 ), 1 )
*
*+    TEST 1
*           Factor A as L*D*L' and compute the ratio
*              norm(L*D*L' - A) / (n * norm(A) * EPS )
*
            CALL SPTTRF( N, D( N+1 ), E( N+1 ), INFO )
*
*           Check error code from SPTTRF.
*
            IF( INFO.NE.IZERO ) THEN
               CALL ALAERH( PATH, 'SPTTRF', INFO, IZERO, ' ', N, N, -1,
     $                      -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 100
            END IF
*
            IF( INFO.GT.0 ) THEN
               RCONDC = ZERO
               GO TO 90
            END IF
*
            CALL SPTT01( N, D, E, D( N+1 ), E( N+1 ), WORK,
     $                   RESULT( 1 ) )
*
*           Print the test ratio if greater than or equal to THRESH.
*
            IF( RESULT( 1 ).GE.THRESH ) THEN
               IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $            CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )N, IMAT, 1, RESULT( 1 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
*
*           Compute RCONDC = 1 / (norm(A) * norm(inv(A))
*
*           Compute norm(A).
*
            ANORM = SLANST( '1', N, D, E )
*
*           Use SPTTRS to solve for one column at a time of inv(A),
*           computing the maximum column sum as we go.
*
            AINVNM = ZERO
            DO 50 I = 1, N
               DO 40 J = 1, N
                  X( J ) = ZERO
   40          CONTINUE
               X( I ) = ONE
               CALL SPTTRS( N, 1, D( N+1 ), E( N+1 ), X, LDA, INFO )
               AINVNM = MAX( AINVNM, SASUM( N, X, 1 ) )
   50       CONTINUE
            RCONDC = ONE / MAX( ONE, ANORM*AINVNM )
*
            DO 80 IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
*
*           Generate NRHS random solution vectors.
*
               IX = 1
               DO 60 J = 1, NRHS
                  CALL SLARNV( 2, ISEED, N, XACT( IX ) )
                  IX = IX + LDA
   60          CONTINUE
*
*           Set the right hand side.
*
               CALL SLAPTM( N, NRHS, ONE, D, E, XACT, LDA, ZERO, B,
     $                      LDA )
*
*+    TEST 2
*           Solve A*x = b and compute the residual.
*
               CALL SLACPY( 'Full', N, NRHS, B, LDA, X, LDA )
               CALL SPTTRS( N, NRHS, D( N+1 ), E( N+1 ), X, LDA, INFO )
*
*           Check error code from SPTTRS.
*
               IF( INFO.NE.0 )
     $            CALL ALAERH( PATH, 'SPTTRS', INFO, 0, ' ', N, N, -1,
     $                         -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
               CALL SLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
               CALL SPTT02( N, NRHS, D, E, X, LDA, WORK, LDA,
     $                      RESULT( 2 ) )
*
*+    TEST 3
*           Check solution from generated exact solution.
*
               CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                      RESULT( 3 ) )
*
*+    TESTS 4, 5, and 6
*           Use iterative refinement to improve the solution.
*
               SRNAMT = 'SPTRFS'
               CALL SPTRFS( N, NRHS, D, E, D( N+1 ), E( N+1 ), B, LDA,
     $                      X, LDA, RWORK, RWORK( NRHS+1 ), WORK, INFO )
*
*           Check error code from SPTRFS.
*
               IF( INFO.NE.0 )
     $            CALL ALAERH( PATH, 'SPTRFS', INFO, 0, ' ', N, N, -1,
     $                         -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
*
               CALL SGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC,
     $                      RESULT( 4 ) )
               CALL SPTT05( N, NRHS, D, E, B, LDA, X, LDA, XACT, LDA,
     $                      RWORK, RWORK( NRHS+1 ), RESULT( 5 ) )
*
*           Print information about the tests that did not pass the
*           threshold.
*
               DO 70 K = 2, 6
                  IF( RESULT( K ).GE.THRESH ) THEN
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                  CALL ALAHD( NOUT, PATH )
                     WRITE( NOUT, FMT = 9998 )N, NRHS, IMAT, K,
     $                  RESULT( K )
                     NFAIL = NFAIL + 1
                  END IF
   70          CONTINUE
               NRUN = NRUN + 5
   80       CONTINUE
*
*+    TEST 7
*           Estimate the reciprocal of the condition number of the
*           matrix.
*
   90       CONTINUE
            SRNAMT = 'SPTCON'
            CALL SPTCON( N, D( N+1 ), E( N+1 ), ANORM, RCOND, RWORK,
     $                   INFO )
*
*           Check error code from SPTCON.
*
            IF( INFO.NE.0 )
     $         CALL ALAERH( PATH, 'SPTCON', INFO, 0, ' ', N, N, -1, -1,
     $                      -1, IMAT, NFAIL, NERRS, NOUT )
*
            RESULT( 7 ) = SGET06( RCOND, RCONDC )
*
*           Print the test ratio if greater than or equal to THRESH.
*
            IF( RESULT( 7 ).GE.THRESH ) THEN
               IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $            CALL ALAHD( NOUT, PATH )
               WRITE( NOUT, FMT = 9999 )N, IMAT, 7, RESULT( 7 )
               NFAIL = NFAIL + 1
            END IF
            NRUN = NRUN + 1
  100    CONTINUE
  110 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( ' N =', I5, ', type ', I2, ', test ', I2, ', ratio = ',
     $      G12.5 )
 9998 FORMAT( ' N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2,
     $      ') = ', G12.5 )
      RETURN
*
*     End of SCHKPT
*
      END
