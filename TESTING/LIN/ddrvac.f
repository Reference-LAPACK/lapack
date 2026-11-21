*> \brief \b DDRVAC
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX,
*                          A, AFAC, B, X, WORK,
*                          RWORK, SWORK, NOUT )
*
*       .. Scalar Arguments ..
*       INTEGER            NMAX, NM, NNS, NOUT
*       DOUBLE PRECISION   THRESH
*       ..
*       .. Array Arguments ..
*       LOGICAL            DOTYPE( * )
*       INTEGER            MVAL( * ), NSVAL( * )
*       REAL               SWORK(*)
*       DOUBLE PRECISION   A( * ), AFAC( * ), B( * ),
*      $                   RWORK( * ), WORK( * ), X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DDRVAC tests DSPOSV.
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
*>          The number of values of N contained in the vector MVAL.
*> \endverbatim
*>
*> \param[in] MVAL
*> \verbatim
*>          MVAL is INTEGER array, dimension (NM)
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
*>          THRESH is DOUBLE PRECISION
*>          The threshold value for the test ratios.  A result is
*>          included in the output file if RESULT >= THRESH.  To have
*>          every test ratio printed, use THRESH = 0.
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
*>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AFAC
*> \verbatim
*>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension
*>                      (NMAX*max(3,NSMAX))
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is DOUBLE PRECISION array, dimension
*>                      (max(2*NMAX,2*NSMAX+NWORK))
*> \endverbatim
*>
*> \param[out] SWORK
*> \verbatim
*>          SWORK is REAL array, dimension
*>                      (NMAX*(NSMAX+NMAX))
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX,
     $                   A, AFAC, B, X, WORK,
     $                   RWORK, SWORK, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            NMAX, NM, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      INTEGER            MVAL( * ), NSVAL( * )
      REAL               SWORK(*)
      DOUBLE PRECISION   A( * ), AFAC( * ), B( * ),
     $                   RWORK( * ), WORK( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 9 )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, TYPE, UPLO, XTYPE
      CHARACTER*3        PATH
      INTEGER            I, IM, IMAT, INFO, IOFF, IRHS, IUPLO,
     $                   IZERO, KL, KU, LDA, MODE, N,
     $                   NERRS, NFAIL, NIMAT, NRHS, NRUN
      DOUBLE PRECISION   ANORM, CNDNUM
*     ..
*     .. Local Arrays ..
      CHARACTER          UPLOS( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. Local Variables ..
      INTEGER            ITER, KASE
      INTEGER            TI, TINFO, TITER
      LOGICAL            TBAD
      DOUBLE PRECISION   TA( 2, 2 ), TB( 2, 1 ), TW( 2, 1 ), TX( 2, 1 )
      REAL               TSW( 6 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, DLACPY, DLARHS, DLATB4, DLATMS,
     $                   DPOT06, DSPOSV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
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
      KASE = 0
      PATH( 1: 1 ) = 'Double precision'
      PATH( 2: 3 ) = 'PO'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      INFOT = 0
*
*     Do for each value of N in MVAL
*
      DO 120 IM = 1, NM
         N = MVAL( IM )
         LDA = MAX( N, 1 )
         NIMAT = NTYPES
         IF( N.LE.0 )
     $      NIMAT = 1
*
         DO 110 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 110
*
*           Skip types 3, 4, or 5 if the matrix size is too small.
*
            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 )
     $         GO TO 110
*
*           Do first for UPLO = 'U', then for UPLO = 'L'
*
            DO 100 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
*
*              Set up parameters with DLATB4 and generate a test matrix
*              with DLATMS.
*
               CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE,
     $                      CNDNUM, DIST )
*
               SRNAMT = 'DLATMS'
               CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE,
     $                      CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK,
     $                      INFO )
*
*              Check error code from DLATMS.
*
               IF( INFO.NE.0 ) THEN
                  CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1,
     $                         -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               END IF
*
*              For types 3-5, zero one row and column of the matrix to
*              test that INFO is returned correctly.
*
               IF( ZEROT ) THEN
                  IF( IMAT.EQ.3 ) THEN
                     IZERO = 1
                  ELSE IF( IMAT.EQ.4 ) THEN
                     IZERO = N
                  ELSE
                     IZERO = N / 2 + 1
                  END IF
                  IOFF = ( IZERO-1 )*LDA
*
*                 Set row and column IZERO of A to 0.
*
                  IF( IUPLO.EQ.1 ) THEN
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
   30                CONTINUE
                  ELSE
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  END IF
               ELSE
                  IZERO = 0
               END IF
*
               DO 60 IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
                  XTYPE = 'N'
*
*                 Form an exact solution and set the right hand side.
*
                  SRNAMT = 'DLARHS'
                  CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU,
     $                         NRHS, A, LDA, X, LDA, B, LDA,
     $                         ISEED, INFO )
*
*                 Compute the L*L' or U'*U factorization of the
*                 matrix and solve the system.
*
                  SRNAMT = 'DSPOSV '
                  KASE = KASE + 1
*
                  CALL DLACPY( 'All', N, N, A, LDA, AFAC, LDA)
*
                  CALL DSPOSV( UPLO, N, NRHS, AFAC, LDA, B, LDA, X, LDA,
     $                         WORK, SWORK, ITER, INFO )

                  IF (ITER.LT.0) THEN
                     CALL DLACPY( 'All', N, N, A, LDA, AFAC, LDA )
                  ENDIF
*
*                 Check error code from DSPOSV .
*
                  IF( INFO.NE.IZERO ) THEN
*
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $                  CALL ALAHD( NOUT, PATH )
                     NERRS = NERRS + 1
*
                     IF( INFO.NE.IZERO .AND. IZERO.NE.0 ) THEN
                        WRITE( NOUT, FMT = 9988 )'DSPOSV',INFO,IZERO,N,
     $                     IMAT
                     ELSE
                        WRITE( NOUT, FMT = 9975 )'DSPOSV',INFO,N,IMAT
                     END IF
                  END IF
*
*                 Skip the remaining test if the matrix is singular.
*
                  IF( INFO.NE.0 )
     $               GO TO 110
*
*                 Check the quality of the solution
*
                  CALL DLACPY( 'All', N, NRHS, B, LDA, WORK, LDA )
*
                  CALL DPOT06( UPLO, N, NRHS, A, LDA, X, LDA, WORK,
     $               LDA, RWORK, RESULT( 1 ) )
*
*                 Check if the test passes the testing.
*                 Print information about the tests that did not
*                 pass the testing.
*
*                 If iterative refinement has been used and claimed to
*                 be successful (ITER>0), we want
*                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1
*
*                 If double precision has been used (ITER<0), we want
*                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
*                 (Cf. the linear solver testing routines)
*
                  IF ((THRESH.LE.0.0E+00)
     $               .OR.((ITER.GE.0).AND.(N.GT.0)
     $               .AND.(RESULT(1).GE.SQRT(DBLE(N))))
     $               .OR.((ITER.LT.0).AND.(RESULT(1).GE.THRESH))) THEN
*
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) THEN
                        WRITE( NOUT, FMT = 8999 )'DPO'
                        WRITE( NOUT, FMT = '( '' Matrix types:'' )' )
                        WRITE( NOUT, FMT = 8979 )
                        WRITE( NOUT, FMT = '( '' Test ratios:'' )' )
                        WRITE( NOUT, FMT = 8960 )1
                        WRITE( NOUT, FMT = '( '' Messages:'' )' )
                     END IF
*
                     WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, 1,
     $                  RESULT( 1 )
*
                     NFAIL = NFAIL + 1
*
                  END IF
*
                  NRUN = NRUN + 1
*
   60          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
*
*     Check that a single precision solve that overflows is
*     detected and that the routine falls back to double precision.
*     For A = diag( 1e-30, 1e-30 ) and B = 1e30, the solution X = 1e60
*     is representable in double but not in single precision, so the
*     single precision solve returns Inf.  We want ITER < 0, INFO = 0
*     and X correct to double precision.
*
      TA( 1, 1 ) = 1.0D-30
      TA( 2, 1 ) = ZERO
      TA( 1, 2 ) = ZERO
      TA( 2, 2 ) = 1.0D-30
      TB( 1, 1 ) = 1.0D+30
      TB( 2, 1 ) = 1.0D+30
      SRNAMT = 'DSPOSV'
      CALL DSPOSV( 'U', 2, 1, TA, 2, TB, 2, TX, 2, TW, TSW,
     $             TITER, TINFO )
      TBAD = TINFO.NE.0 .OR. TITER.GE.0
      DO 130 TI = 1, 2
         IF( .NOT.( ABS( TX( TI, 1 )-1.0D+60 ).LE.1.0D+50 ) )
     $      TBAD = .TRUE.
  130 CONTINUE
      NRUN = NRUN + 1
      IF( TBAD ) THEN
         IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 )
     $      CALL ALAHD( NOUT, PATH )
         WRITE( NOUT, FMT = 9987 )'DSPOSV', TITER, TINFO
         NFAIL = NFAIL + 1
      END IF
*
*     Print a summary of the results.
*
      IF( NFAIL.GT.0 ) THEN
         WRITE( NOUT, FMT = 9996 )'DSPOSV', NFAIL, NRUN
      ELSE
         WRITE( NOUT, FMT = 9995 )'DSPOSV', NRUN
      END IF
      IF( NERRS.GT.0 ) THEN
         WRITE( NOUT, FMT = 9994 )NERRS
      END IF
*
 9998 FORMAT( ' UPLO=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9987 FORMAT( ' *** ', A6, ' did not fall back to double precision',
     $      ' when the single precision solve overflowed: ITER =',
     $      I5, ', INFO =', I5 )
 9996 FORMAT( 1X, A6, ': ', I6, ' out of ', I6,
     $      ' tests failed to pass the threshold' )
 9995 FORMAT( /1X, 'All tests for ', A6,
     $      ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( 6X, I6, ' error messages recorded' )
*
*     SUBNAM, INFO, INFOE, N, IMAT
*
 9988 FORMAT( ' *** ', A6, ' returned with INFO =', I5, ' instead of ',
     $      I5, / ' ==> N =', I5, ', type ',
     $      I2 )
*
*     SUBNAM, INFO, N, IMAT
*
 9975 FORMAT( ' *** Error code from ', A6, '=', I5, ' for M=', I5,
     $      ', type ', I2 )
 8999 FORMAT( / 1X, A3, ':  positive definite dense matrices' )
 8979 FORMAT( 4X, '1. Diagonal', 24X, '7. Last n/2 columns zero', / 4X,
     $      '2. Upper triangular', 16X,
     $      '8. Random, CNDNUM = sqrt(0.1/EPS)', / 4X,
     $      '3. Lower triangular', 16X, '9. Random, CNDNUM = 0.1/EPS',
     $      / 4X, '4. Random, CNDNUM = 2', 13X,
     $      '10. Scaled near underflow', / 4X, '5. First column zero',
     $      14X, '11. Scaled near overflow', / 4X,
     $      '6. Last column zero' )
 8960 FORMAT( 3X, I2, ': norm_1( B - A * X )  / ',
     $      '( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF',
     $      / 4x, 'or norm_1( B - A * X )  / ',
     $      '( norm_1(A) * norm_1(X) * EPS ) > THRES if DPOTRF' )

      RETURN
*
*     End of DDRVAC
*
      END
