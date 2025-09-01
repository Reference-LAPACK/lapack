*> \brief \b CGET35
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE CGET35( RMAX, LMAX, NINFO, KNT, NIN )
*
*       .. Scalar Arguments ..
*       INTEGER            KNT, LMAX, NIN, NINFO
*       REAL               RMAX
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CGET35 tests CTRSYL, a routine for solving the Sylvester matrix
*> equation
*>
*>    op(A)*X + ISGN*X*op(B) = scale*C,
*>
*> A and B are assumed to be in Schur canonical form, op() represents an
*> optional transpose, and ISGN can be -1 or +1.  Scale is an output
*> less than or equal to 1, chosen to avoid overflow in X.
*>
*> The test code verifies that the following residual is order 1:
*>
*>    norm(op(A)*X + ISGN*X*op(B) - scale*C) /
*>        (EPS*max(norm(A),norm(B))*norm(X))
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[out] RMAX
*> \verbatim
*>          RMAX is REAL
*>          Value of the largest test ratio.
*> \endverbatim
*>
*> \param[out] LMAX
*> \verbatim
*>          LMAX is INTEGER
*>          Example number where largest test ratio achieved.
*> \endverbatim
*>
*> \param[out] NINFO
*> \verbatim
*>          NINFO is INTEGER
*>          Number of examples where INFO is nonzero.
*> \endverbatim
*>
*> \param[out] KNT
*> \verbatim
*>          KNT is INTEGER
*>          Total number of examples tested.
*> \endverbatim
*>
*> \param[in] NIN
*> \verbatim
*>          NIN is INTEGER
*>          Input logical unit number.
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
      SUBROUTINE CGET35( RMAX, LMAX, NINFO, KNT, NIN )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            KNT, LMAX, NIN, NINFO
      REAL               RMAX
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            LDT
      PARAMETER          ( LDT = 10 )
      REAL               ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0E0, ONE = 1.0E0, TWO = 2.0E0 )
      REAL               LARGE
      PARAMETER          ( LARGE = 1.0E6 )
      COMPLEX            CONE
      PARAMETER          ( CONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANA, TRANB
      INTEGER            I, IMLA, IMLAD, IMLB, IMLC, INFO, ISGN, ITRANA,
     $                   ITRANB, J, M, N
      REAL               BIGNUM, EPS, RES, RES1, SCALE, SMLNUM, TNRM,
     $                   XNRM
      COMPLEX            RMUL
*     ..
*     .. Local Arrays ..
      REAL               DUM( 1 ), VM1( 3 ), VM2( 3 )
      COMPLEX            A( LDT, LDT ), ATMP( LDT, LDT ), B( LDT, LDT ),
     $                   BTMP( LDT, LDT ), C( LDT, LDT ),
     $                   CSAV( LDT, LDT ), CTMP( LDT, LDT )
*     ..
*     .. External Functions ..
      REAL               CLANGE, SLAMCH
      EXTERNAL           CLANGE, SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CTRSYL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, REAL, SQRT
*     ..
*     .. Executable Statements ..
*
*     Get machine parameters
*
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Set up test case parameters
*
      VM1( 1 ) = SQRT( SMLNUM )
      VM1( 2 ) = ONE
      VM1( 3 ) = LARGE
      VM2( 1 ) = ONE
      VM2( 2 ) = ONE + TWO*EPS
      VM2( 3 ) = TWO
*
      KNT = 0
      NINFO = 0
      LMAX = 0
      RMAX = ZERO
*
*     Begin test loop
*
   10 CONTINUE
      READ( NIN, FMT = * )M, N
      IF( N.EQ.0 )
     $   RETURN
      DO 20 I = 1, M
         READ( NIN, FMT = * )( ATMP( I, J ), J = 1, M )
   20 CONTINUE
      DO 30 I = 1, N
         READ( NIN, FMT = * )( BTMP( I, J ), J = 1, N )
   30 CONTINUE
      DO 40 I = 1, M
         READ( NIN, FMT = * )( CTMP( I, J ), J = 1, N )
   40 CONTINUE
      DO 170 IMLA = 1, 3
         DO 160 IMLAD = 1, 3
            DO 150 IMLB = 1, 3
               DO 140 IMLC = 1, 3
                  DO 130 ITRANA = 1, 2
                     DO 120 ITRANB = 1, 2
                        DO 110 ISGN = -1, 1, 2
                           IF( ITRANA.EQ.1 )
     $                        TRANA = 'N'
                           IF( ITRANA.EQ.2 )
     $                        TRANA = 'C'
                           IF( ITRANB.EQ.1 )
     $                        TRANB = 'N'
                           IF( ITRANB.EQ.2 )
     $                        TRANB = 'C'
                           TNRM = ZERO
                           DO 60 I = 1, M
                              DO 50 J = 1, M
                                 A( I, J ) = ATMP( I, J )*VM1( IMLA )
                                 TNRM = MAX( TNRM, ABS( A( I, J ) ) )
   50                         CONTINUE
                              A( I, I ) = A( I, I )*VM2( IMLAD )
                              TNRM = MAX( TNRM, ABS( A( I, I ) ) )
   60                      CONTINUE
                           DO 80 I = 1, N
                              DO 70 J = 1, N
                                 B( I, J ) = BTMP( I, J )*VM1( IMLB )
                                 TNRM = MAX( TNRM, ABS( B( I, J ) ) )
   70                         CONTINUE
   80                      CONTINUE
                           IF( TNRM.EQ.ZERO )
     $                        TNRM = ONE
                           DO 100 I = 1, M
                              DO 90 J = 1, N
                                 C( I, J ) = CTMP( I, J )*VM1( IMLC )
                                 CSAV( I, J ) = C( I, J )
   90                         CONTINUE
  100                      CONTINUE
                           KNT = KNT + 1
                           CALL CTRSYL( TRANA, TRANB, ISGN, M, N, A,
     $                                  LDT, B, LDT, C, LDT, SCALE,
     $                                  INFO )
                           IF( INFO.NE.0 )
     $                        NINFO = NINFO + 1
                           XNRM = CLANGE( 'M', M, N, C, LDT, DUM )
                           RMUL = CONE
                           IF( XNRM.GT.ONE .AND. TNRM.GT.ONE ) THEN
                              IF( XNRM.GT.BIGNUM / TNRM ) THEN
                                 RMUL = MAX( XNRM, TNRM )
                                 RMUL = CONE / RMUL
                              END IF
                           END IF
                           CALL CGEMM( TRANA, 'N', M, N, M, RMUL, A,
     $                                 LDT, C, LDT, -SCALE*RMUL, CSAV,
     $                                 LDT )
                           CALL CGEMM( 'N', TRANB, M, N, N,
     $                                 REAL( ISGN )*RMUL, C, LDT, B,
     $                                 LDT, CONE, CSAV, LDT )
                           RES1 = CLANGE( 'M', M, N, CSAV, LDT, DUM )
                           RES = RES1 / MAX( SMLNUM, SMLNUM*XNRM,
     $                           ( ( ABS( RMUL )*TNRM )*EPS )*XNRM )
                           IF( RES.GT.RMAX ) THEN
                              LMAX = KNT
                              RMAX = RES
                           END IF
  110                   CONTINUE
  120                CONTINUE
  130             CONTINUE
  140          CONTINUE
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE
      GO TO 10
*
*     End of CGET35
*
      END
