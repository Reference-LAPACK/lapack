*> \brief \b SGBEQU
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SGBEQU + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgbequ.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgbequ.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgbequ.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND,
*                          AMAX, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, KL, KU, LDAB, M, N
*       REAL               AMAX, COLCND, ROWCND
*       ..
*       .. Array Arguments ..
*       REAL               AB( LDAB, * ), C( * ), R( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SGBEQU computes row and column scalings intended to equilibrate an
*> M-by-N band matrix A and reduce its condition number.  R returns the
*> row scale factors and C the column scale factors, chosen to try to
*> make the largest element in each row and column of the matrix B with
*> elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.
*>
*> R(i) and C(j) are restricted to be between SMLNUM = smallest safe
*> number and BIGNUM = largest safe number.  Use of these scaling
*> factors is not guaranteed to reduce the condition number of A but
*> works well in practice.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] KL
*> \verbatim
*>          KL is INTEGER
*>          The number of subdiagonals within the band of A.  KL >= 0.
*> \endverbatim
*>
*> \param[in] KU
*> \verbatim
*>          KU is INTEGER
*>          The number of superdiagonals within the band of A.  KU >= 0.
*> \endverbatim
*>
*> \param[in] AB
*> \verbatim
*>          AB is REAL array, dimension (LDAB,N)
*>          The band matrix A, stored in rows 1 to KL+KU+1.  The j-th
*>          column of A is stored in the j-th column of the array AB as
*>          follows:
*>          AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl).
*> \endverbatim
*>
*> \param[in] LDAB
*> \verbatim
*>          LDAB is INTEGER
*>          The leading dimension of the array AB.  LDAB >= KL+KU+1.
*> \endverbatim
*>
*> \param[out] R
*> \verbatim
*>          R is REAL array, dimension (M)
*>          If INFO = 0, or INFO > M, R contains the row scale factors
*>          for A.
*> \endverbatim
*>
*> \param[out] C
*> \verbatim
*>          C is REAL array, dimension (N)
*>          If INFO = 0, C contains the column scale factors for A.
*> \endverbatim
*>
*> \param[out] ROWCND
*> \verbatim
*>          ROWCND is REAL
*>          If INFO = 0 or INFO > M, ROWCND contains the ratio of the
*>          smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and
*>          AMAX is neither too large nor too small, it is not worth
*>          scaling by R.
*> \endverbatim
*>
*> \param[out] COLCND
*> \verbatim
*>          COLCND is REAL
*>          If INFO = 0, COLCND contains the ratio of the smallest
*>          C(i) to the largest C(i).  If COLCND >= 0.1, it is not
*>          worth scaling by C.
*> \endverbatim
*>
*> \param[out] AMAX
*> \verbatim
*>          AMAX is REAL
*>          Absolute value of largest matrix element.  If AMAX is very
*>          close to overflow or very close to underflow, the matrix
*>          should be scaled.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*>          > 0:  if INFO = i, and i is
*>                <= M:  the i-th row of A is exactly zero
*>                >  M:  the (i-M)-th column of A is exactly zero
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
*> \ingroup gbequ
*
*  =====================================================================
      SUBROUTINE SGBEQU( M, N, KL, KU, AB, LDAB, R, C, ROWCND,
     $                   COLCND,
     $                   AMAX, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KL, KU, LDAB, M, N
      REAL               AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      REAL               AB( LDAB, * ), C( * ), R( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, KD
      REAL               BIGNUM, RCMAX, RCMIN, SMLNUM
*     ..
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KL+KU+1 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGBEQU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      END IF
*
*     Get machine constants.
*
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
*     Compute row scale factors.
*
      DO 10 I = 1, M
         R( I ) = ZERO
   10 CONTINUE
*
*     Find the maximum element in each row.
*
      KD = KU + 1
      DO 30 J = 1, N
         DO 20 I = MAX( J-KU, 1 ), MIN( J+KL, M )
            R( I ) = MAX( R( I ), ABS( AB( KD+I-J, J ) ) )
   20    CONTINUE
   30 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 50 I = 1, M
            IF( R( I ).EQ.ZERO ) THEN
               INFO = I
               RETURN
            END IF
   50    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE
*
*        Compute ROWCND = min(R(I)) / max(R(I))
*
         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
*     Compute column scale factors
*
      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE
*
*     Find the maximum element in each column,
*     assuming the row scaling computed above.
*
      KD = KU + 1
      DO 90 J = 1, N
         DO 80 I = MAX( J-KU, 1 ), MIN( J+KL, M )
            C( J ) = MAX( C( J ), ABS( AB( KD+I-J, J ) )*R( I ) )
   80    CONTINUE
   90 CONTINUE
*
*     Find the maximum and minimum scale factors.
*
      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE
*
      IF( RCMIN.EQ.ZERO ) THEN
*
*        Find the first zero scale factor and return an error code.
*
         DO 110 J = 1, N
            IF( C( J ).EQ.ZERO ) THEN
               INFO = M + J
               RETURN
            END IF
  110    CONTINUE
      ELSE
*
*        Invert the scale factors.
*
         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE
*
*        Compute COLCND = min(C(J)) / max(C(J))
*
         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      END IF
*
      RETURN
*
*     End of SGBEQU
*
      END
