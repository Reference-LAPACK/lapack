*> \brief \b SORT01
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE SORT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RESID )
*
*       .. Scalar Arguments ..
*       CHARACTER          ROWCOL
*       INTEGER            LDU, LWORK, M, N
*       REAL               RESID
*       ..
*       .. Array Arguments ..
*       REAL               U( LDU, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SORT01 checks that the matrix U is orthogonal by computing the ratio
*>
*>    RESID = norm( I - U*U' ) / ( n * EPS ), if ROWCOL = 'R',
*> or
*>    RESID = norm( I - U'*U ) / ( m * EPS ), if ROWCOL = 'C'.
*>
*> Alternatively, if there isn't sufficient workspace to form
*> I - U*U' or I - U'*U, the ratio is computed as
*>
*>    RESID = abs( I - U*U' ) / ( n * EPS ), if ROWCOL = 'R',
*> or
*>    RESID = abs( I - U'*U ) / ( m * EPS ), if ROWCOL = 'C'.
*>
*> where EPS is the machine precision.  ROWCOL is used only if m = n;
*> if m > n, ROWCOL is assumed to be 'C', and if m < n, ROWCOL is
*> assumed to be 'R'.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ROWCOL
*> \verbatim
*>          ROWCOL is CHARACTER
*>          Specifies whether the rows or columns of U should be checked
*>          for orthogonality.  Used only if M = N.
*>          = 'R':  Check for orthogonal rows of U
*>          = 'C':  Check for orthogonal columns of U
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix U.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix U.
*> \endverbatim
*>
*> \param[in] U
*> \verbatim
*>          U is REAL array, dimension (LDU,N)
*>          The orthogonal matrix U.  U is checked for orthogonal columns
*>          if m > n or if m = n and ROWCOL = 'C'.  U is checked for
*>          orthogonal rows if m < n or if m = n and ROWCOL = 'R'.
*> \endverbatim
*>
*> \param[in] LDU
*> \verbatim
*>          LDU is INTEGER
*>          The leading dimension of the array U.  LDU >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (LWORK)
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The length of the array WORK.  For best performance, LWORK
*>          should be at least N*(N+1) if ROWCOL = 'C' or M*(M+1) if
*>          ROWCOL = 'R', but the test will be done even if LWORK is 0.
*> \endverbatim
*>
*> \param[out] RESID
*> \verbatim
*>          RESID is REAL
*>          RESID = norm( I - U * U' ) / ( n * EPS ), if ROWCOL = 'R', or
*>          RESID = norm( I - U' * U ) / ( m * EPS ), if ROWCOL = 'C'.
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
*> \ingroup single_eig
*
*  =====================================================================
      SUBROUTINE SORT01( ROWCOL, M, N, U, LDU, WORK, LWORK, RESID )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ROWCOL
      INTEGER            LDU, LWORK, M, N
      REAL               RESID
*     ..
*     .. Array Arguments ..
      REAL               U( LDU, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANSU
      INTEGER            I, J, K, LDWORK, MNMIN
      REAL               EPS, TMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT, SLAMCH, SLANSY
      EXTERNAL           LSAME, SDOT, SLAMCH, SLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASET, SSYRK
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      RESID = ZERO
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      EPS = SLAMCH( 'Precision' )
      IF( M.LT.N .OR. ( M.EQ.N .AND. LSAME( ROWCOL, 'R' ) ) ) THEN
         TRANSU = 'N'
         K = N
      ELSE
         TRANSU = 'T'
         K = M
      END IF
      MNMIN = MIN( M, N )
*
      IF( ( MNMIN+1 )*MNMIN.LE.LWORK ) THEN
         LDWORK = MNMIN
      ELSE
         LDWORK = 0
      END IF
      IF( LDWORK.GT.0 ) THEN
*
*        Compute I - U*U' or I - U'*U.
*
         CALL SLASET( 'Upper', MNMIN, MNMIN, ZERO, ONE, WORK, LDWORK )
         CALL SSYRK( 'Upper', TRANSU, MNMIN, K, -ONE, U, LDU, ONE, WORK,
     $               LDWORK )
*
*        Compute norm( I - U*U' ) / ( K * EPS ) .
*
         RESID = SLANSY( '1', 'Upper', MNMIN, WORK, LDWORK,
     $           WORK( LDWORK*MNMIN+1 ) )
         RESID = ( RESID / REAL( K ) ) / EPS
      ELSE IF( TRANSU.EQ.'T' ) THEN
*
*        Find the maximum element in abs( I - U'*U ) / ( m * EPS )
*
         DO 20 J = 1, N
            DO 10 I = 1, J
               IF( I.NE.J ) THEN
                  TMP = ZERO
               ELSE
                  TMP = ONE
               END IF
               TMP = TMP - SDOT( M, U( 1, I ), 1, U( 1, J ), 1 )
               RESID = MAX( RESID, ABS( TMP ) )
   10       CONTINUE
   20    CONTINUE
         RESID = ( RESID / REAL( M ) ) / EPS
      ELSE
*
*        Find the maximum element in abs( I - U*U' ) / ( n * EPS )
*
         DO 40 J = 1, M
            DO 30 I = 1, J
               IF( I.NE.J ) THEN
                  TMP = ZERO
               ELSE
                  TMP = ONE
               END IF
               TMP = TMP - SDOT( N, U( J, 1 ), LDU, U( I, 1 ), LDU )
               RESID = MAX( RESID, ABS( TMP ) )
   30       CONTINUE
   40    CONTINUE
         RESID = ( RESID / REAL( N ) ) / EPS
      END IF
      RETURN
*
*     End of SORT01
*
      END
