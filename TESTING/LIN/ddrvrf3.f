*> \brief \b DDRVRF3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DDRVRF3( NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2,
*      +                    D_WORK_DLANGE, D_WORK_DGEQRF, TAU )
*
*       .. Scalar Arguments ..
*       INTEGER            LDA, NN, NOUT
*       DOUBLE PRECISION   THRESH
*       ..
*       .. Array Arguments ..
*       INTEGER            NVAL( NN )
*       DOUBLE PRECISION   A( LDA, * ), ARF( * ), B1( LDA, * ),
*      +                   B2( LDA, * ), D_WORK_DGEQRF( * ),
*      +                   D_WORK_DLANGE( * ), TAU( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DDRVRF3 tests the LAPACK RFP routines:
*>     DTFSM
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>                The unit number for output.
*> \endverbatim
*>
*> \param[in] NN
*> \verbatim
*>          NN is INTEGER
*>                The number of values of N contained in the vector NVAL.
*> \endverbatim
*>
*> \param[in] NVAL
*> \verbatim
*>          NVAL is INTEGER array, dimension (NN)
*>                The values of the matrix dimension N.
*> \endverbatim
*>
*> \param[in] THRESH
*> \verbatim
*>          THRESH is DOUBLE PRECISION
*>                The threshold value for the test ratios.  A result is
*>                included in the output file if RESULT >= THRESH.  To have
*>                every test ratio printed, use THRESH = 0.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,NMAX)
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>                The leading dimension of the array A.  LDA >= max(1,NMAX).
*> \endverbatim
*>
*> \param[out] ARF
*> \verbatim
*>          ARF is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2).
*> \endverbatim
*>
*> \param[out] B1
*> \verbatim
*>          B1 is DOUBLE PRECISION array, dimension (LDA,NMAX)
*> \endverbatim
*>
*> \param[out] B2
*> \verbatim
*>          B2 is DOUBLE PRECISION array, dimension (LDA,NMAX)
*> \endverbatim
*>
*> \param[out] D_WORK_DLANGE
*> \verbatim
*>          D_WORK_DLANGE is DOUBLE PRECISION array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] D_WORK_DGEQRF
*> \verbatim
*>          D_WORK_DGEQRF is DOUBLE PRECISION array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (NMAX)
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
      SUBROUTINE DDRVRF3( NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2,
     +                    D_WORK_DLANGE, D_WORK_DGEQRF, TAU )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            LDA, NN, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            NVAL( NN )
      DOUBLE PRECISION   A( LDA, * ), ARF( * ), B1( LDA, * ),
     +                   B2( LDA, * ), D_WORK_DGEQRF( * ),
     +                   D_WORK_DLANGE( * ), TAU( * )
*     ..
*
*  =====================================================================
*     ..
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) ,
     +                     ONE  = ( 1.0D+0, 0.0D+0 ) )
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 1 )
*     ..
*     .. Local Scalars ..
      CHARACTER          UPLO, CFORM, DIAG, TRANS, SIDE
      INTEGER            I, IFORM, IIM, IIN, INFO, IUPLO, J, M, N, NA,
     +                   NFAIL, NRUN, ISIDE, IDIAG, IALPHA, ITRANS
      DOUBLE PRECISION   EPS, ALPHA
*     ..
*     .. Local Arrays ..
      CHARACTER          UPLOS( 2 ), FORMS( 2 ), TRANSS( 2 ),
     +                   DIAGS( 2 ), SIDES( 2 )
      INTEGER            ISEED( 4 ), ISEEDY( 4 )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE, DLARND
      EXTERNAL           DLAMCH, DLANGE, DLARND, LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRTTF, DGEQRF, DGEQLF, DTFSM, DTRSM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Scalars in Common ..
      CHARACTER*32       SRNAMT
*     ..
*     .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS  / 'U', 'L' /
      DATA               FORMS  / 'N', 'T' /
      DATA               SIDES  / 'L', 'R' /
      DATA               TRANSS / 'N', 'T' /
      DATA               DIAGS  / 'N', 'U' /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      NRUN = 0
      NFAIL = 0
      INFO = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = DLAMCH( 'Precision' )
*
      DO 170 IIM = 1, NN
*
         M = NVAL( IIM )
*
         DO 160 IIN = 1, NN
*
            N = NVAL( IIN )
*
            DO 150 IFORM = 1, 2
*
               CFORM = FORMS( IFORM )
*
               DO 140 IUPLO = 1, 2
*
                  UPLO = UPLOS( IUPLO )
*
                  DO 130 ISIDE = 1, 2
*
                     SIDE = SIDES( ISIDE )
*
                     DO 120 ITRANS = 1, 2
*
                        TRANS = TRANSS( ITRANS )
*
                        DO 110 IDIAG = 1, 2
*
                           DIAG = DIAGS( IDIAG )
*
                           DO 100 IALPHA = 1, 3
*
                              IF ( IALPHA.EQ.1 ) THEN
                                 ALPHA = ZERO
                              ELSE IF ( IALPHA.EQ.2 ) THEN
                                 ALPHA = ONE
                              ELSE
                                 ALPHA = DLARND( 2, ISEED )
                              END IF
*
*                             All the parameters are set:
*                                CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
*                                and ALPHA
*                             READY TO TEST!
*
                              NRUN = NRUN + 1
*
                              IF ( ISIDE.EQ.1 ) THEN
*
*                                The case ISIDE.EQ.1 is when SIDE.EQ.'L'
*                                -> A is M-by-M ( B is M-by-N )
*
                                 NA = M
*
                              ELSE
*
*                                The case ISIDE.EQ.2 is when SIDE.EQ.'R'
*                                -> A is N-by-N ( B is M-by-N )
*
                                 NA = N
*
                              END IF
*
*                             Generate A our NA--by--NA triangular
*                             matrix.
*                             Our test is based on forward error so we
*                             do want A to be well conditioned! To get
*                             a well-conditioned triangular matrix, we
*                             take the R factor of the QR/LQ factorization
*                             of a random matrix.
*
                              DO J = 1, NA
                                 DO I = 1, NA
                                    A( I, J ) = DLARND( 2, ISEED )
                                 END DO
                              END DO
*
                              IF ( IUPLO.EQ.1 ) THEN
*
*                                The case IUPLO.EQ.1 is when SIDE.EQ.'U'
*                                -> QR factorization.
*
                                 SRNAMT = 'DGEQRF'
                                 CALL DGEQRF( NA, NA, A, LDA, TAU,
     +                                        D_WORK_DGEQRF, LDA,
     +                                        INFO )
*
*                                Forcing main diagonal of test matrix to
*                                be unit makes it ill-conditioned for
*                                some test cases
*
                                 IF ( LSAME( DIAG, 'U' ) ) THEN
                                    DO J = 1, NA
                                       DO I = 1, J
                                          A( I, J ) = A( I, J ) /
     +                                            ( 2.0 * A( J, J ) )
                                       END DO
                                    END DO
                                 END IF
*
                              ELSE
*
*                                The case IUPLO.EQ.2 is when SIDE.EQ.'L'
*                                -> QL factorization.
*
                                 SRNAMT = 'DGELQF'
                                 CALL DGELQF( NA, NA, A, LDA, TAU,
     +                                        D_WORK_DGEQRF, LDA,
     +                                        INFO )
*
*                                Forcing main diagonal of test matrix to
*                                be unit makes it ill-conditioned for
*                                some test cases
*
                                 IF ( LSAME( DIAG, 'U' ) ) THEN
                                    DO I = 1, NA
                                       DO J = 1, I
                                          A( I, J ) = A( I, J ) /
     +                                            ( 2.0 * A( I, I ) )
                                       END DO
                                    END DO
                                 END IF
*
                              END IF
*
*                             Store a copy of A in RFP format (in ARF).
*
                              SRNAMT = 'DTRTTF'
                              CALL DTRTTF( CFORM, UPLO, NA, A, LDA, ARF,
     +                                     INFO )
*
*                             Generate B1 our M--by--N right-hand side
*                             and store a copy in B2.
*
                              DO J = 1, N
                                 DO I = 1, M
                                    B1( I, J ) = DLARND( 2, ISEED )
                                    B2( I, J ) = B1( I, J )
                                 END DO
                              END DO
*
*                             Solve op( A ) X = B or X op( A ) = B
*                             with DTRSM
*
                              SRNAMT = 'DTRSM'
                              CALL DTRSM( SIDE, UPLO, TRANS, DIAG, M, N,
     +                               ALPHA, A, LDA, B1, LDA )
*
*                             Solve op( A ) X = B or X op( A ) = B
*                             with DTFSM
*
                              SRNAMT = 'DTFSM'
                              CALL DTFSM( CFORM, SIDE, UPLO, TRANS,
     +                                    DIAG, M, N, ALPHA, ARF, B2,
     +                                    LDA )
*
*                             Check that the result agrees.
*
                              DO J = 1, N
                                 DO I = 1, M
                                    B1( I, J ) = B2( I, J ) - B1( I, J )
                                 END DO
                              END DO
*
                              RESULT( 1 ) = DLANGE( 'I', M, N, B1, LDA,
     +                                            D_WORK_DLANGE )
*
                              RESULT( 1 ) = RESULT( 1 ) / SQRT( EPS )
     +                                    / MAX ( MAX( M, N ), 1 )
*
                              IF( RESULT( 1 ).GE.THRESH ) THEN
                                 IF( NFAIL.EQ.0 ) THEN
                                    WRITE( NOUT, * )
                                    WRITE( NOUT, FMT = 9999 )
                                 END IF
                                 WRITE( NOUT, FMT = 9997 ) 'DTFSM',
     +                              CFORM, SIDE, UPLO, TRANS, DIAG, M,
     +                              N, RESULT( 1 )
                                 NFAIL = NFAIL + 1
                              END IF
*
  100                      CONTINUE
  110                   CONTINUE
  120                CONTINUE
  130             CONTINUE
  140          CONTINUE
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE
*
*     Print a summary of the results.
*
      IF ( NFAIL.EQ.0 ) THEN
         WRITE( NOUT, FMT = 9996 ) 'DTFSM', NRUN
      ELSE
         WRITE( NOUT, FMT = 9995 ) 'DTFSM', NFAIL, NRUN
      END IF
*
 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing DTFSM
     +         ***')
 9997 FORMAT( 1X, '     Failure in ',A5,', CFORM=''',A1,''',',
     + ' SIDE=''',A1,''',',' UPLO=''',A1,''',',' TRANS=''',A1,''',',
     + ' DIAG=''',A1,''',',' M=',I3,', N =', I3,', test=',G12.5)
 9996 FORMAT( 1X, 'All tests for ',A5,' auxiliary routine passed the ',
     +        'threshold ( ',I5,' tests run)')
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5,
     +        ' tests failed to pass the threshold')
*
      RETURN
*
*     End of DDRVRF3
*
      END
