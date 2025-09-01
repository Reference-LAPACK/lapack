*> \brief \b CCKGSV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
*                          NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R,
*                          IWORK, WORK, RWORK, NIN, NOUT, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, NIN, NM, NMATS, NMAX, NOUT
*       REAL               THRESH
*       ..
*       .. Array Arguments ..
*       INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ),
*      $                   PVAL( * )
*       REAL               ALPHA( * ), BETA( * ), RWORK( * )
*       COMPLEX            A( * ), AF( * ), B( * ), BF( * ), Q( * ),
*      $                   R( * ), U( * ), V( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CCKGSV tests CGGSVD:
*>        the GSVD for M-by-N matrix A and P-by-N matrix B.
*> \endverbatim
*
*  Arguments:
*  ==========
*
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
*> \param[in] PVAL
*> \verbatim
*>          PVAL is INTEGER array, dimension (NP)
*>          The values of the matrix row dimension P.
*> \endverbatim
*>
*> \param[in] NVAL
*> \verbatim
*>          NVAL is INTEGER array, dimension (NN)
*>          The values of the matrix column dimension N.
*> \endverbatim
*>
*> \param[in] NMATS
*> \verbatim
*>          NMATS is INTEGER
*>          The number of matrix types to be tested for each combination
*>          of matrix dimensions.  If NMATS >= NTYPES (the maximum
*>          number of matrix types), then all the different types are
*>          generated for testing.  If NMATS < NTYPES, another input line
*>          is read to get the numbers of the matrix types to be used.
*> \endverbatim
*>
*> \param[in,out] ISEED
*> \verbatim
*>          ISEED is INTEGER array, dimension (4)
*>          On entry, the seed of the random number generator.  The array
*>          elements should be between 0 and 4095, otherwise they will be
*>          reduced mod 4096, and ISEED(4) must be odd.
*>          On exit, the next seed in the random number sequence after
*>          all the test matrices have been generated.
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
*> \param[in] NMAX
*> \verbatim
*>          NMAX is INTEGER
*>          The maximum value permitted for M or N, used in dimensioning
*>          the work arrays.
*> \endverbatim
*>
*> \param[out] A
*> \verbatim
*>          A is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] AF
*> \verbatim
*>          AF is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] B
*> \verbatim
*>          B is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] BF
*> \verbatim
*>          BF is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] U
*> \verbatim
*>          U is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] V
*> \verbatim
*>          V is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] ALPHA
*> \verbatim
*>          ALPHA is REAL array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] BETA
*> \verbatim
*>          BETA is REAL array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (NMAX)
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (NMAX*NMAX)
*> \endverbatim
*>
*> \param[out] RWORK
*> \verbatim
*>          RWORK is REAL array, dimension (NMAX)
*> \endverbatim
*>
*> \param[in] NIN
*> \verbatim
*>          NIN is INTEGER
*>          The unit number for input.
*> \endverbatim
*>
*> \param[in] NOUT
*> \verbatim
*>          NOUT is INTEGER
*>          The unit number for output.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0 :  successful exit
*>          > 0 :  If CLATMS returns an error code, the absolute value
*>                 of it is returned.
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
*> \ingroup complex_eig
*
*  =====================================================================
      SUBROUTINE CCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH,
     $                   NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R,
     $                   IWORK, WORK, RWORK, NIN, NOUT, INFO )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, NIN, NM, NMATS, NMAX, NOUT
      REAL               THRESH
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ),
     $                   PVAL( * )
      REAL               ALPHA( * ), BETA( * ), RWORK( * )
      COMPLEX            A( * ), AF( * ), B( * ), BF( * ), Q( * ),
     $                   R( * ), U( * ), V( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NTESTS
      PARAMETER          ( NTESTS = 12 )
      INTEGER            NTYPES
      PARAMETER          ( NTYPES = 8 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRSTT
      CHARACTER          DISTA, DISTB, TYPE
      CHARACTER*3        PATH
      INTEGER            I, IINFO, IM, IMAT, KLA, KLB, KUA, KUB, LDA,
     $                   LDB, LDQ, LDR, LDU, LDV, LWORK, M, MODEA,
     $                   MODEB, N, NFAIL, NRUN, NT, P, K, L
      REAL               ANORM, BNORM, CNDNMA, CNDNMB
*     ..
*     .. Local Arrays ..
      LOGICAL            DOTYPE( NTYPES )
      REAL               RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAHDG, ALAREQ, ALASUM, CLATMS, SLATB9, CGSVTS3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      PATH( 1: 3 ) = 'GSV'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDA = NMAX
      LDB = NMAX
      LDU = NMAX
      LDV = NMAX
      LDQ = NMAX
      LDR = NMAX
      LWORK = NMAX*NMAX
*
*     Specific cases
*
*     Test: https://github.com/Reference-LAPACK/lapack/issues/411#issue-608776973
*
      M = 6
      P = 6
      N = 6
      A(1:M*N) = CMPLX(1.E0, 0.E0)
      B(1:M*N) = CMPLX(0.E0, 0.E0)
      B(1+0*M) = CMPLX(9.E19, 0.E0)
      B(2+1*M) = CMPLX(9.E18, 0.E0)
      B(3+2*M) = CMPLX(9.E17, 0.E0)
      B(4+3*M) = CMPLX(9.E16, 0.E0)
      B(5+4*M) = CMPLX(9.E15, 0.E0)
      B(6+5*M) = CMPLX(9.E14, 0.E0)
      CALL CGGSVD3('N','N','N', M, P, N, K, L, A, M, B, M,
     $              ALPHA, BETA, U, 1, V, 1, Q, 1,
     $              WORK, M*N, RWORK, IWORK, INFO)
*
*     Print information there is a NAN in BETA
      DO 40 I = 1, L
         IF( BETA(I).NE.BETA(I) ) THEN
            INFO = -I
            EXIT
         END IF
   40 CONTINUE
      IF( INFO.LT.0 ) THEN
         IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
            FIRSTT = .FALSE.
            CALL ALAHDG( NOUT, PATH )
         END IF
         WRITE( NOUT, FMT = 9997 ) -INFO
         NFAIL = NFAIL + 1
      END IF
      NRUN = NRUN + 1
      INFO = 0
*
*     Do for each value of M in MVAL.
*
      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         N = NVAL( IM )
*
         DO 20 IMAT = 1, NTYPES
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) )
     $         GO TO 20
*
*           Set up parameters with SLATB9 and generate test
*           matrices A and B with CLATMS.
*
            CALL SLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB,
     $                   ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB,
     $                   DISTA, DISTB )
*
*           Generate M by N matrix A
*
            CALL CLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA,
     $                   ANORM, KLA, KUA, 'No packing', A, LDA, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF
*
*           Generate P by N matrix B
*
            CALL CLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB,
     $                   BNORM, KLB, KUB, 'No packing', B, LDB, WORK,
     $                   IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF
*
            NT = 6
*
            CALL CGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V,
     $                    LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK,
     $                    LWORK, RWORK, RESULT )
*
*           Print information about the tests that did not
*           pass the threshold.
*
            DO 10 I = 1, NT
               IF( RESULT( I ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                     FIRSTT = .FALSE.
                     CALL ALAHDG( NOUT, PATH )
                  END IF
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I,
     $               RESULT( I )
                  NFAIL = NFAIL + 1
               END IF
   10       CONTINUE
            NRUN = NRUN + NT
*
   20    CONTINUE
   30 CONTINUE
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )
*
 9999 FORMAT( ' CLATMS in CCKGSV   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' FOUND NaN in BETA(', I4,')' )
      RETURN
*
*     End of CCKGSV
*
      END
