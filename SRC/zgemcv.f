*> \brief \b ZGEMCV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZGEMCV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*       .. Scalar Arguments ..
*       COMPLEX*16 ALPHA,BETA
*       INTEGER INCX,INCY,LDA,M,N
*       CHARACTER TRANS
*       ..
*       .. Array Arguments ..
*       COMPLEX*16 A(LDA,*),X(*),Y(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZGEMCV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*conjg(x) + beta*y,   or   y := alpha*A**T*conjg(x) + beta*y,   or
*>
*>    y := alpha*A**H*conjg(x) + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   y := alpha*A*conjg(x) + beta*y.
*>
*>              TRANS = 'T' or 't'   y := alpha*A**T*conjg(x) + beta*y.
*>
*>              TRANS = 'C' or 'c'   y := alpha*A**H*conjg(x) + beta*y.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension ( LDA, N )
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is COMPLEX*16
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is COMPLEX*16 array, dimension at least
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*>           If either m or n is zero, then Y not referenced and the function
*>           performs a quick return.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
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
*> \ingroup gemv
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZGEMCV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
      IMPLICIT NONE
*
*       .. Scalar Arguments ..
        COMPLEX*16   ALPHA,BETA
        INTEGER      INCX,INCY,LDA,M,N
        CHARACTER    TRANS
*       ..
*       .. Array Arguments ..
        COMPLEX*16   A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
        COMPLEX*16   ZERO
        PARAMETER    (ZERO = (0.0D+0, 0.0D+0))
*     ..
*     .. Local Scalars ..
        LOGICAL      NOTRANS, NOCONJ
        INTEGER      KX, KY, I, J, IY, JX, INFO, LENX, LENY
        COMPLEX*16   TEMP
*     ..
*     .. External Functions ..
        LOGICAL      LSAME
        EXTERNAL     LSAME
*     ..
*     .. External Subroutines ..
        EXTERNAL     ZSCAL
*     ..
*     .. Intrinsic Functions ..
        INTRINSIC    MAX, CONJG
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.(LSAME(TRANS, 'N').OR.LSAME(TRANS, 'T').OR.
     $   LSAME(TRANS,'C'))) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('ZGEMCV',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF (M.EQ.0.OR.N.EQ.0.OR.(ALPHA.EQ.ZERO.AND.BETA.EQ.ZERO)) THEN
         RETURN
      END IF
      NOTRANS = LSAME(TRANS, 'N')
      NOCONJ  = LSAME(TRANS, 'T')
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (NOTRANS) THEN
*
*         Note: this means A is LENY by LENX
*
          LENX = N
          LENY = M
      ELSE
*
*         Note: this means A is LENX by LENY
*
          LENX = M
          LENY = N
      END IF
*
*     Note: I do not use the increments less than 0 cases so these
*     are not tested, however the general access pattern follows that
*     of zgemv, so it should work
*
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Compute y = beta*y
*     zscal should do nothing if beta is 1, so no extra check is needed
*     at this level
*
      CALL ZSCAL(LENY, BETA, Y, INCY)
*
*     If alpha is 0, we are done
*
      IF( ALPHA.EQ.ZERO ) THEN
         RETURN
      END IF
*
*     Now, we compute y = alpha*op(A)*x + y one component of y at a time
*     our formulation tries to let the compiler take advantage of FMAs
*     and takes advantage of the fact that if op(A) = A, A is LENY by LENX
*     and otherwise, A is LENX by LENY to simplify the loops
*
*     We do not currently have different loops for the differing cases
*     of INCX and INCY. This can be changed if desired.
*
      IF( NOTRANS ) THEN
*
*        op(A) = A
*
         IY = KY
         DO I = 1, LENY
            TEMP = ZERO
            JX = KX
            DO J = 1, LENX
               TEMP = TEMP + A(I,J) * CONJG(X(JX))
               JX = JX + INCX
            END DO
            Y(IY) = Y(IY) + TEMP*ALPHA
            IY = IY + INCY
         END DO
      ELSE IF(NOCONJ) THEN
*
*        op(A) = A**T
*
         IY = KY
         DO I = 1, LENY
            TEMP = ZERO
            JX = KX
            DO J = 1, LENX
               TEMP = TEMP + A(J,I) * CONJG(X(JX))
               JX = JX + INCX
            END DO
            Y(IY) = Y(IY) + TEMP*ALPHA
            IY = IY + INCY
         END DO
      ELSE
*
*        op(A) = A**H
*
         IY = KY
         DO I = 1, LENY
            TEMP = ZERO
            JX = KX
            DO J = 1, LENX
               TEMP = TEMP + CONJG(A(J,I) * X(JX))
               JX = JX + INCX
            END DO
            Y(IY) = Y(IY) + TEMP*ALPHA
            IY = IY + INCY
         END DO
      END IF
      END SUBROUTINE
