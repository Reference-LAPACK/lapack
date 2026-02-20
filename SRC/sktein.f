*> \brief \b SKTEIN
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SKTEIN + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sktein.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sktein.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sktein.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKTEIN( N, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
*                          IWORK, IFAIL, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, LDZ, M, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
*      $                   IWORK( * )
*       REAL               E( * ), W( * ), WORK( * ), Z( LDZ, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKTEIN computes the eigenvectors of a real skew-symmetric tridiagonal
*> matrix T corresponding to specified eigenvalues, using inverse
*> iteration.
*>
*> The maximum number of iterations allowed for each eigenvector is
*> specified by an internal parameter MAXITS (currently set to 5).
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix.  N >= 0.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is REAL array, dimension (N-1)
*>          The (n-1) subdiagonal elements of the tridiagonal matrix
*>          T, in elements 1 to N-1.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of eigenvectors to be found.  0 <= M <= N.
*> \endverbatim
*>
*> \param[in] W
*> \verbatim
*>          W is REAL array, dimension (N)
*>          The first M elements of W contain the eigenvalues for
*>          which eigenvectors are to be computed.  The eigenvalues
*>          should be grouped by split-off block and ordered from
*>          smallest to largest within the block.  ( The output array
*>          W from SKTEBZ with ORDER = 'B' is expected here. )
*> \endverbatim
*>
*> \param[in] IBLOCK
*> \verbatim
*>          IBLOCK is INTEGER array, dimension (N)
*>          The submatrix indices associated with the corresponding
*>          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
*>          the first submatrix from the top, =2 if W(i) belongs to
*>          the second submatrix, etc.  ( The output array IBLOCK
*>          from SSTEBZ is expected here. )
*> \endverbatim
*>
*> \param[in] ISPLIT
*> \verbatim
*>          ISPLIT is INTEGER array, dimension (N)
*>          The splitting points, at which T breaks up into submatrices.
*>          The first submatrix consists of rows/columns 1 to
*>          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
*>          through ISPLIT( 2 ), etc.
*>          ( The output array ISPLIT from SKTEBZ is expected here. )
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ, M)
*>          The computed eigenvectors.  The eigenvector associated
*>          with the eigenvalue W(i) is stored in the i-th column of
*>          Z.  Any vector which fails to converge is set to its current
*>          iterate after MAXITS iterations.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.  LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (6*N)
*> \endverbatim
*>
*> \param[out] IWORK
*> \verbatim
*>          IWORK is INTEGER array, dimension (N+1)
*> \endverbatim
*>
*> \param[out] IFAIL
*> \verbatim
*>          IFAIL is INTEGER array, dimension (M)
*>          On normal exit, all elements of IFAIL are zero.
*>          If one or more eigenvectors fail to converge after
*>          MAXITS iterations, then their indices are stored in
*>          array IFAIL.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit.
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, then i eigenvectors failed to converge
*>               in MAXITS iterations.  Their indices are stored in
*>               array IFAIL.
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  MAXITS  INTEGER, default = 5
*>          The maximum number of iterations performed.
*>
*>  EXTRA   INTEGER, default = 2
*>          The number of iterations performed after norm growth
*>          criterion is satisfied, should be at least 1.
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
*> \ingroup stein
*
*  =====================================================================
      SUBROUTINE SKTEIN( N, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,
     $                   IWORK, IFAIL, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDZ, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IBLOCK( * ), IFAIL( * ), ISPLIT( * ),
     $                   IWORK( * )
      REAL               E( * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE, TWO, FOUR, TEN, ODM3, ODM1
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0,
     $                   TWO = 2.0E+0, FOUR = 4.0E+0, TEN = 1.0E+1,
     $                   ODM3 = 1.0E-3, ODM1 = 1.0E-1 )
      INTEGER            MAXITS, EXTRA
      PARAMETER          ( MAXITS = 5, EXTRA = 2 )
*     ..
*     .. Local Scalars ..
      INTEGER            B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1,
     $                   INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1,
     $                   JBLK, JMAX, JMAX1, JMAX2, NBLK, NRMCHK
      REAL               CTR, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL,
     $                   SCL, SEP, STPCRT, TOL, XJ, XJM, S1, S2, K1,
     $                   K2, P1, P2, P3, Q1, Q2, Q3, TS1, TS2, TK1,
     $                   TK2
*     ..
*     .. Local Arrays ..
      INTEGER            ISEED( 8 )
*     ..
*     .. External Functions ..
      INTEGER            ISAMAX
      REAL               SDOT, SLAMCH, SNRM2
      EXTERNAL           ISAMAX, SDOT, SLAMCH, SNRM2
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SLAGTF, SLAGTS, SLAGTFK,
     $                   SLAGTSK, SLARNV, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      DO 10 I = 1, M
         IFAIL( I ) = 0
   10 CONTINUE
*
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( M.LT.0 .OR. M.GT.N ) THEN
         INFO = -3
      ELSE IF( LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE
         DO 20 J = 1, M
            IF( J.GT.1 .AND. IBLOCK( J ).LT.IBLOCK( J-1 ) ) THEN
               INFO = -5
               GO TO 30
            END IF
            IF( J.GT.2 .AND. W( J ).NE.ZERO .AND.
     $          W( J-2 ).NE.ZERO .AND.
     $          IBLOCK( J ).EQ.IBLOCK( J-2 ) .AND.
     $          W( J ).GT.W( J-2 ) ) THEN
               INFO = -4
               GO TO 30
            END IF
            IF ( W( J ).LT.ZERO ) THEN
               INFO = -4
               GO TO 30
            END IF
   20    CONTINUE
   30    CONTINUE
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKTEIN', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 ) THEN
         RETURN
      ELSE IF( N.EQ.1 ) THEN
         Z( 1, 1 ) = ONE
         
      ELSE IF( N.EQ.2 ) THEN
         IF( E(1).GE.ZERO ) THEN
            Z( 1, 1 ) = ONE
            Z( 1, 2 ) = ZERO
            Z( 2, 1 ) = ZERO
            Z( 2, 2 ) = ONE
         ELSE
            Z( 1, 1 ) = ONE
            Z( 1, 2 ) = ZERO
            Z( 2, 1 ) = ZERO
            Z( 2, 2 ) = -ONE
         END IF
         RETURN
      END IF
*
*     Get machine constants.
*
      EPS = SLAMCH( 'Precision' )
*
*     Initialize seed for random number generator SLARNV.
*
      DO 40 I = 1, 8
         ISEED( I ) = 1
   40 CONTINUE
*
*     Initialize pointers.
*
      INDRV1 = 0
      INDRV2 = INDRV1 + 2 * N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
*
*     Compute eigenvectors of matrix blocks.
*
      J1 = 1
      DO 230 NBLK = 1, IBLOCK( M )
*
*        Find starting and ending indices of block nblk.
*
         IF( NBLK.EQ.1 ) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT( NBLK-1 ) + 1
         END IF
         BN = ISPLIT( NBLK )
         BLKSIZ = BN - B1 + 1
         IF( BLKSIZ.EQ.1 .OR. BLKSIZ.EQ.2 )
     $      GO TO 60
         GPIND = J1
*
*        Compute reorthogonalization criterion and stopping criterion.
*
         ONENRM = MAX(ABS( E( B1 ) ), ABS( E( BN-1 ) ) )
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX( ONENRM, ABS( E( I-1 ) ) +
     $                    ABS( E( I ) ) )
   50    CONTINUE
         ORTOL = ODM3*ONENRM
*
         STPCRT = SQRT( ODM1 / REAL( BLKSIZ ) )
*
*        Loop through eigenvalues of block nblk.
*
   60    CONTINUE
         JBLK = 0
         J = J1
         DO WHILE (J.LE.M)
            IF( IBLOCK( J ).NE.NBLK ) THEN
               J1 = J
               GO TO 230
            END IF
*
            IF (W( J ).NE.ZERO) THEN
*
*              A pair of conjugate eigenvalues.
*
               JBLK = JBLK + 2
               XJ = W( J )
*
*              Skip all the work if the block size is one.
*
               IF( BLKSIZ.EQ.2 ) THEN
                  WORK( INDRV1+1 ) = ONE
                  WORK( INDRV1+2 ) = ZERO
                  WORK( N+INDRV1+1 ) = ZERO
                  WORK( N+INDRV1+2 ) = ONE
                  GO TO 120
               END IF
*
*              If eigenvalues j and j-1 are too close, add a relatively
*              small perturbation.
*
               IF( JBLK.GT.2 ) THEN
                  EPS1 = ABS( EPS*XJ )
                  PERTOL = TEN*EPS1
                  SEP = XJM - XJ
                  IF( SEP.LT.PERTOL )
     $               XJ = XJM - PERTOL
               END IF
*
               ITS = 0
               NRMCHK = 0
*
*              Get random starting vector.
*
               CALL SLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
               CALL SLARNV( 2, ISEED+4, BLKSIZ,
     $                      WORK( N+INDRV1+1 ) )
*
*              Copy the matrix T so it won't be destroyed in factorization.
*
               DO I = 1, BLKSIZ
                  WORK( INDRV4+I ) = ZERO
               END DO
               CALL SCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ),
     $                     1 )
               CALL SSCAL( BLKSIZ-1, -ONE, WORK( INDRV2+2 ), 1 )
               CALL SCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ),
     $                     1 )
*
*              Compute LU factors with partial pivoting  ( PT = LU )
*
               TOL = ZERO
               CALL SLAGTFK( BLKSIZ, WORK( INDRV4+1 ), -XJ,
     $                       WORK( INDRV2+2 ), WORK( INDRV3+1 ),
     $                       TOL, WORK( INDRV5+1 ), IWORK, IINFO )
*
*              Update iteration count.
*
   70          CONTINUE
               ITS = ITS + 1
               IF( ITS.GT.MAXITS )
     $            GO TO 100
*
*              Normalize and scale the righthand side vector Pb.
*
               JMAX1 = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               JMAX2 = ISAMAX( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
               SCL = REAL( BLKSIZ )*ONENRM*MAX( EPS,
     $               ABS( WORK( INDRV4+BLKSIZ ) ) ) /
     $               MAX( ABS( WORK( INDRV1+JMAX1 ) ),
     $               ABS( WORK( INDRV1+JMAX2 ) ))
               CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
               CALL SSCAL( BLKSIZ, SCL, WORK( N+INDRV1+1 ), 1 )
*
*              Solve the system LU = Pb.
*
               CALL SLAGTSK( -1, BLKSIZ, WORK( INDRV4+1 ),
     $                       WORK( INDRV2+2 ), WORK( INDRV3+1 ),
     $                       WORK( INDRV5+1 ), IWORK,
     $                       WORK( INDRV1+1 ), WORK( N+INDRV1+1 ),
     $                       TOL, IINFO )
*
*              Equalize the L2 norm of vectors in a pair
*
               JMAX1 = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               JMAX2 = ISAMAX( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
               S1 = SNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
               S2 = SNRM2( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
               IF( S1.GT.S2 ) THEN
                  SCL = S2 / S1
                  CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
               ELSE
                  SCL = S1 / S2
                  CALL SSCAL( BLKSIZ, SCL, WORK( N+INDRV1+1 ), 1 )
               END IF
               IF ( ABS( WORK( INDRV1+JMAX1 ) ).GT.
     $              ABS( WORK( N+INDRV1+JMAX2 ) ) .AND.
     $              WORK( INDRV1+JMAX1 ).LT.ZERO .OR.
     $              ABS( WORK( INDRV1+JMAX1 ) ).LE.
     $              ABS( WORK( N+INDRV1+JMAX2 ) ) .AND.
     $              WORK( N+INDRV1+JMAX2 ).LT.ZERO ) THEN
                  CALL SSCAL( BLKSIZ, -ONE, WORK( INDRV1+1 ), 1 )
                  CALL SSCAL( BLKSIZ, -ONE, WORK( N+INDRV1+1 ), 1 )
               END IF
               CALL SCOPY( BLKSIZ, WORK( INDRV1+1 ), 1,
     $                     Z( B1, J ), 1 )
               CALL SCOPY( BLKSIZ, WORK( N+INDRV1+1 ), 1,
     $                     Z( B1, J+1 ), 1 )
               CALL SAXPY( BLKSIZ, -ONE, Z( B1, J+1 ), 1,
     $                     WORK( INDRV1+1 ), 1 )
               CALL SAXPY( BLKSIZ, ONE, Z( B1, J ), 1,
     $                     WORK( N+INDRV1+1 ), 1 )
               SCL = SQRT( TWO ) / TWO
               CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
               CALL SSCAL( BLKSIZ, SCL, WORK( N+INDRV1+1 ), 1 )
               IF ( XJ.LT.ORTOL / TWO ) THEN
                  S1 = SNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
                  S2 = SNRM2( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
                  IF( S1.GT.S2 ) THEN
                     SCL = S2 / S1
                     CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
                  ELSE
                     SCL = S1 / S2
                     CALL SSCAL( BLKSIZ, SCL, WORK( N+INDRV1+1 ), 1 )
                  END IF
               END IF
*
*              Reorthogonalize by modified Gram-Schmidt if eigenvalues are
*              close enough.
*
               IF( JBLK.EQ.2 )
     $            GO TO 90
               IF( ABS( XJM-XJ ).GT.ORTOL )
     $            GPIND = J
               IF( GPIND.NE.J ) THEN
                  DO 80 I = GPIND, J - 2, 2
                     CALL SCOPY( BLKSIZ, WORK( INDRV1+1 ), 1,
     $                           Z( B1, J ), 1 )
                     CALL SCOPY( BLKSIZ, WORK( N+INDRV1+1 ), 1,
     $                           Z( B1, J+1 ), 1 )
                     S1 = SDOT( BLKSIZ, WORK( INDRV1+1 ), 1,
     $                          Z( B1, I ), 1 )
                     K1 = SDOT( BLKSIZ, WORK( N+INDRV1+1 ), 1,
     $                          Z( B1, I ), 1 )
                     S2 = SDOT( BLKSIZ, WORK( N+INDRV1+1 ), 1,
     $                          Z( B1, I+1 ), 1 )
                     K2 = SDOT( BLKSIZ, WORK( INDRV1+1 ), 1,
     $                          Z( B1, I+1 ), 1 )
                     SCL = MAX( MAX(ABS(S1), ABS(K1)), MAX(ABS(S2),
     $                         ABS(K2)) )
                     TS1 = S1 / SCL;
                     TK1 = K1 / SCL;
                     TS2 = S2 / SCL;
                     TK2 = K2 / SCL;
                     IF( ABS(TS1*TK1 + TS2*TK2) .GT. ZERO ) THEN
                        CTR = (TK1*TK1 - TK2*TK2 + TS2*TS2 -
     $                         TS1*TS1) / (TS1*TK1 + TS2*TK2)
                        IF( CTR.LT.ZERO ) THEN
                           P1 = (SQRT(CTR*CTR + FOUR) - CTR) / TWO
                        ELSE
                           P1 = (-SQRT(CTR*CTR + FOUR) - CTR) / TWO
                        END IF
                        Q1 = -P1
                        P2 = S1 - P1 * K1
                        Q2 = S2 - Q1 * K2
                        P3 = K2 - P1 * S2
                        Q3 = K1 - Q1 * S1
                        CALL SAXPY( BLKSIZ, -P1, Z( B1, J+1 ), 1,
     $                              WORK( INDRV1+1 ), 1 )
                        CALL SAXPY( BLKSIZ, -Q1, Z( B1, J ), 1,
     $                              WORK( N+INDRV1+1 ), 1 )
                        CALL SAXPY( BLKSIZ, -P2, Z( B1, I ), 1,
     $                              WORK( INDRV1+1 ), 1 )
                        CALL SAXPY( BLKSIZ, -Q2, Z( B1, I+1 ), 1,
     $                              WORK( N+INDRV1+1 ), 1 )
                        CALL SAXPY( BLKSIZ, -P3, Z( B1, I+1 ), 1,
     $                              WORK( INDRV1+1 ), 1 )
                        CALL SAXPY( BLKSIZ, -Q3, Z( B1, I ), 1,
     $                              WORK( N+INDRV1+1 ), 1 )
                        IF ( XJ.LT.ORTOL / TWO ) THEN
                           S1 = SNRM2( BLKSIZ, WORK( INDRV1+1 ),
     $                                 1 )
                           S2 = SNRM2( BLKSIZ, WORK( N+INDRV1+1 ),
     $                                 1 )
                           IF( S1.GT.S2 ) THEN
                              SCL = S2 / S1
                              CALL SSCAL( BLKSIZ, SCL,
     $                                    WORK( INDRV1+1 ), 1 )
                           ELSE
                              SCL = S1 / S2
                              CALL SSCAL( BLKSIZ, SCL,
     $                                    WORK( N+INDRV1+1 ), 1 )
                           END IF
                        END IF
                     END IF
   80             CONTINUE
               END IF
*
*              Check the infinity norm of the iterate.
*
   90          CONTINUE
               JMAX1 = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               JMAX2 = ISAMAX( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
               NRM = MAX( ABS( WORK( INDRV1+JMAX1 ) ),
     $               ABS( WORK( N+INDRV1+JMAX2 ) ) )
*
*              Continue for additional iterations after norm reaches
*              stopping criterion.
*
               IF( NRM.LT.STPCRT )
     $            GO TO 70
               NRMCHK = NRMCHK + 1
               IF( NRMCHK.LT.EXTRA+1 )
     $            GO TO 70
*
               GO TO 110
*
*              If stopping criterion was not satisfied, update info and
*              store eigenvector number in array ifail.
*
  100          CONTINUE
               INFO = INFO + 2
               IFAIL( INFO ) = J
*
*              Accept iterate as jth eigenvector.
*
  110          CONTINUE
               JMAX1 = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               SCL = ONE / SNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
               CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
               JMAX2 = ISAMAX( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
               SCL = ONE / SNRM2( BLKSIZ, WORK( N+INDRV1+1 ), 1 )
               CALL SSCAL( BLKSIZ, SCL, WORK( N+INDRV1+1 ), 1 )
               IF ( ABS( WORK( INDRV1+JMAX1 ) ).GT.
     $              ABS( WORK( N+INDRV1+JMAX2 ) ) .AND.
     $              WORK( INDRV1+JMAX1 ).LT.ZERO .OR.
     $              ABS( WORK( INDRV1+JMAX1 ) ).LE.
     $              ABS( WORK( N+INDRV1+JMAX2 ) ) .AND.
     $              WORK( N+INDRV1+JMAX2 ).LT.ZERO ) THEN
                  CALL SSCAL( BLKSIZ, -ONE, WORK( INDRV1+1 ), 1 )
                  CALL SSCAL( BLKSIZ, -ONE, WORK( N+INDRV1+1 ), 1 )
               END IF
  120          CONTINUE
               DO 130 I = 1, N
                  Z( I, J ) = ZERO
                  Z( I, J+1 ) = ZERO
  130          CONTINUE
               DO 140 I = 1, BLKSIZ
                  Z( B1+I-1, J ) = WORK( INDRV1+I )
                  Z( B1+I-1, J+1 ) = WORK( N+INDRV1+I )
  140          CONTINUE
*
*              Save the shift to check eigenvalue spacing at next
*              iteration.
*
               XJM = XJ
*
               J = J + 2
*
            ELSE
*
*              Single eigenvalues (zero).
*
               JBLK = JBLK + 1
               XJ = W( J )
*
*              Skip all the work if the block size is one.
*
               IF( BLKSIZ.EQ.1 ) THEN
                  WORK( INDRV1+1 ) = ONE
                  GO TO 200
               END IF
*
*              If eigenvalues j and j-1 are too close, add a relatively
*              small perturbation.
*
               IF( JBLK.GT.1 ) THEN
                  EPS1 = ABS( EPS*XJ )
                  PERTOL = TEN*EPS1
                  SEP = XJM - XJ
                  IF( SEP.LT.PERTOL )
     $               XJ = XJM - PERTOL
               END IF
*
               ITS = 0
               NRMCHK = 0
*
*              Get random starting vector.
*
               CALL SLARNV( 2, ISEED, BLKSIZ, WORK( INDRV1+1 ) )
*
*              Copy the matrix T so it won't be destroyed in factorization.
*
               DO I = 1, BLKSIZ
                  WORK( INDRV4+I ) = ZERO
               END DO
               CALL SCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV2+2 ),
     $                     1 )
               CALL SSCAL( BLKSIZ-1, -ONE, WORK( INDRV2+2 ), 1 )
               CALL SCOPY( BLKSIZ-1, E( B1 ), 1, WORK( INDRV3+1 ),
     $                     1 )
*
*              Compute LU factors with partial pivoting  ( PT = LU )
*
               TOL = ZERO
               CALL SLAGTF( BLKSIZ, WORK( INDRV4+1 ), XJ,
     $                      WORK( INDRV2+2 ), WORK( INDRV3+1 ),
     $                      TOL, WORK( INDRV5+1 ), IWORK,
     $                      IINFO )
*
*              Update iteration count.
*
  150          CONTINUE
               ITS = ITS + 1
               IF( ITS.GT.MAXITS )
     $            GO TO 180
*
*              Normalize and scale the righthand side vector Pb.
*
               JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               SCL = REAL( BLKSIZ )*ONENRM*MAX( EPS,
     $               ABS( WORK( INDRV4+BLKSIZ ) ) ) /
     $               ABS( WORK( INDRV1+JMAX ) )
               CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
*
*              Solve the system LU = Pb.
*
               CALL SLAGTS( -1, BLKSIZ, WORK( INDRV4+1 ),
     $                      WORK( INDRV2+2 ),
     $                      WORK( INDRV3+1 ), WORK( INDRV5+1 ),
     $                      IWORK, WORK( INDRV1+1 ), TOL, IINFO )
*
*              Reorthogonalize by modified Gram-Schmidt if eigenvalues are
*              close enough.
*
               IF( JBLK.EQ.1 )
     $            GO TO 170
               IF( ABS( XJM-XJ ).GT.ORTOL )
     $            GPIND = J
               IF( GPIND.NE.J ) THEN
                  DO 160 I = GPIND, J - 1
                     CTR = -SDOT( BLKSIZ, WORK( INDRV1+1 ), 1,
     $                            Z( B1, I ), 1 )
                     CALL SAXPY( BLKSIZ, CTR, Z( B1, I ), 1,
     $                           WORK( INDRV1+1 ), 1 )
  160             CONTINUE
               END IF
*
*              Check the infinity norm of the iterate.
*
  170          CONTINUE
               JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               NRM = ABS( WORK( INDRV1+JMAX ) )
*
*              Continue for additional iterations after norm reaches
*              stopping criterion.
*
               IF( NRM.LT.STPCRT )
     $            GO TO 150
               NRMCHK = NRMCHK + 1
               IF( NRMCHK.LT.EXTRA+1 )
     $            GO TO 150
*
               GO TO 190
*
*              If stopping criterion was not satisfied, update info and
*              store eigenvector number in array ifail.
*
  180          CONTINUE
               INFO = INFO + 1
               IFAIL( INFO ) = J
*
*              Accept iterate as jth eigenvector.
*
  190          CONTINUE
               SCL = ONE / SNRM2( BLKSIZ, WORK( INDRV1+1 ), 1 )
               JMAX = ISAMAX( BLKSIZ, WORK( INDRV1+1 ), 1 )
               IF( WORK( INDRV1+JMAX ).LT.ZERO )
     $            SCL = -SCL
               CALL SSCAL( BLKSIZ, SCL, WORK( INDRV1+1 ), 1 )
  200          CONTINUE
               DO 210 I = 1, N
                  Z( I, J ) = ZERO
  210          CONTINUE
               DO 220 I = 1, BLKSIZ
                  Z( B1+I-1, J ) = WORK( INDRV1+I )
  220          CONTINUE
*
*              Save the shift to check eigenvalue spacing at next
*              iteration.
*
               XJM = XJ
*
               J = J + 1
            END IF
*
         END DO
  230 CONTINUE
*
      RETURN
*
*     End of SKTEIN
*
      END
