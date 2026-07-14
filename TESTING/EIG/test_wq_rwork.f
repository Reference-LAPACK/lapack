*> \brief Test workspace query precision for z* RWORK
*
*  =========== DOCUMENTATION ===========
*
*  Purpose
*  =======
*
*  TEST_WQ_RWORK validates that workspace query calls (LWORK=-1,
*  LRWORK=-1, LIWORK=-1) return exact RWORK sizes for COMPLEX*16
*  routines with O(N^2) LRWMIN formulas.
*
*  When LRWMIN > 2^24 (approx N > 2896 for formula 1+5N+2N^2),
*  storing the value through a REAL (float32) intermediary loses
*  precision. This test catches that regression by checking
*  INT(RWORK(1)) == expected at N values above the threshold.
*
*  No large matrices are allocated -- workspace queries return
*  immediately after storing sizes, so the test runs in microseconds.
*
*  =========== END DOCUMENTATION ========
*
      PROGRAM TEST_WQ_RWORK
*
      IMPLICIT NONE
*
*     .. Parameters ..
      INTEGER            NNVALS
      PARAMETER          ( NNVALS = 3 )
*
*     .. Local Scalars ..
      INTEGER            INFO, N, LRWEXP, NFAIL, NPASS, I, LRWGOT
*
*     .. Local Arrays ..
*     Minimal dummy arrays for workspace queries (never accessed
*     by the routines when LWORK=-1).
      INTEGER            NVALS( NNVALS ), IWORK( 1 )
      COMPLEX*16         A( 1 ), B( 1 ), AB( 1 ), BB( 1 )
      COMPLEX*16         AP( 1 ), BP( 1 ), Z( 1 ), WORK( 1 )
      DOUBLE PRECISION   W( 1 ), RWORK( 1 ), D( 1 ), E( 1 )
*
*     .. External Subroutines ..
      EXTERNAL           ZHEEVD, ZHEGVD, ZHBEVD, ZHPEVD
      EXTERNAL           ZHPGVD, ZHBGVD, ZSTEDC
*
*     Test N values: 1000 (below 2^24 threshold), 3000 and 5000 (above)
      DATA NVALS / 1000, 3000, 5000 /
*
*     .. Executable Statements ..
*
      NFAIL = 0
      NPASS = 0
*
      WRITE( *, * ) 'Workspace query precision test for z* RWORK'
      WRITE( *, * ) '============================================'
      WRITE( *, * )
*
      DO 100 I = 1, NNVALS
         N = NVALS( I )
*
*        Expected LRWMIN for JOBZ='V': 1 + 5*N + 2*N**2
*        (common to ZHEEVD, ZHEGVD, ZHBEVD, ZHPEVD, ZHPGVD, ZHBGVD)
*
         LRWEXP = 1 + 5*N + 2*N*N
*
*        ---- ZHEEVD ----
*
         INFO = 0
         CALL ZHEEVD( 'V', 'U', N, A, N, W,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZHEEVD ', N, LRWEXP, LRWGOT, INFO
         END IF
*
*        ---- ZHEGVD ----
*
         INFO = 0
         CALL ZHEGVD( 1, 'V', 'U', N, A, N, B, N, W,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZHEGVD ', N, LRWEXP, LRWGOT, INFO
         END IF
*
*        ---- ZHBEVD ----
*        KD=0 (diagonal band matrix), LDAB=1, LDZ=N
*
         INFO = 0
         CALL ZHBEVD( 'V', 'U', N, 0, AB, 1, W, Z, N,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZHBEVD ', N, LRWEXP, LRWGOT, INFO
         END IF
*
*        ---- ZHPEVD ----
*        LDZ=N
*
         INFO = 0
         CALL ZHPEVD( 'V', 'U', N, AP, W, Z, N,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZHPEVD ', N, LRWEXP, LRWGOT, INFO
         END IF
*
*        ---- ZHPGVD ----
*        ITYPE=1, LDZ=N
*
         INFO = 0
         CALL ZHPGVD( 1, 'V', 'U', N, AP, BP, W, Z, N,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZHPGVD ', N, LRWEXP, LRWGOT, INFO
         END IF
*
*        ---- ZHBGVD ----
*        KA=0, KB=0, LDAB=1, LDBB=1, LDZ=N
*
         INFO = 0
         CALL ZHBGVD( 'V', 'U', N, 0, 0, AB, 1, BB, 1,
     $                W, Z, N,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZHBGVD ', N, LRWEXP, LRWGOT, INFO
         END IF
*
*        ---- ZSTEDC (COMPZ='I') ----
*        Expected: 1 + 4*N + 2*N**2
*
         LRWEXP = 1 + 4*N + 2*N*N
*
         INFO = 0
         CALL ZSTEDC( 'I', N, D, E, Z, N,
     $                WORK, -1, RWORK, -1, IWORK, -1, INFO )
         LRWGOT = INT( RWORK( 1 ) )
         IF( INFO.EQ.0 .AND. LRWGOT.EQ.LRWEXP ) THEN
            NPASS = NPASS + 1
         ELSE
            NFAIL = NFAIL + 1
            WRITE( *, 9999 ) 'ZSTEDC ', N, LRWEXP, LRWGOT, INFO
         END IF
*
  100 CONTINUE
*
*     Print summary
*
      WRITE( *, * )
      IF( NFAIL.EQ.0 ) THEN
         WRITE( *, 9998 ) NPASS
      ELSE
         WRITE( *, 9997 ) NFAIL, NFAIL + NPASS
      END IF
*
 9999 FORMAT( ' FAIL: ', A7, ' N=', I6,
     $        ' expected=', I12, ' got=', I12, ' INFO=', I4 )
 9998 FORMAT( ' All ', I3, ' workspace query tests PASSED' )
 9997 FORMAT( ' ', I3, ' of ', I3, ' tests FAILED' )
*
      IF( NFAIL.NE.0 ) STOP 1
*
      END
