*> \brief \b DORGLK computes the explicit Q factor from DGELQF and DLARFT
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE DORGLK(M, N, Q, LDQ)
*
*        .. Scalar Arguments ..
*        INTEGER           M, N, LDQ
*        ..
*        .. Array Arguments ..
*        DOUBLE PRECISION  Q(LDQ,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGLK generates an m by n real matrix Q with orthonormal columns,
*> which is defined as the first n rows of the product of n
*> elementary reflectors
*>
*>       Q  =  I - V'*T*V = H(1) H(2) . . . H(n)
*>
*> Where V is an m by n matrix whose rows are householder reflectors
*> as returned by DGELQF and T is the n by n matrix returned by DLARFT
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix V, and the order of T. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix V. N >= 0.
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>       Q is DOUBLE PRECISION array, dimension (LDQ,N)
*>       On entry, the lower triangular part and diagonal contains
*>       The array T as returned from DLARFT. In addition, the
*>       strictly upper triangular portion of the i-th row contains
*>       the vector which defines the elementary reflector H(i),
*>       for i = 1,2,...,m, as returned by DGELQF
*>       On exit, the m-by-n matrix Q.
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
*  =====================================================================
      SUBROUTINE DORGLK(M, N, Q, LDQ)
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, LDQ
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION  Q(LDQ,*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION  NEG_ONE, ONE
      PARAMETER(NEG_ONE=-1.0D+0, ONE=1.0D+0)
*     ..
*     .. Local Scalars ..
      INTEGER           I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL          DTRMM, DTRTRM, DLUMM
*     ..
*     .. Intrinsic Functions..
      INTRINSIC         MIN
*     ..
*     .. Executable Statements ..
*
*     Break Q apart as follows
*
*           |-----|
*     Q =   | T V |
*           |-----|
*
*     Where T is an m-by-m lower triangular matrix, and V is an
*     m-by-n assumed unit upper trapezoidal matrix
*
*     In turn, break apart V as follows
*
*           |---------|
*     V =   | V_1 V_2 |
*           |---------|
*
*     Where:
*
*     V_1 \in \R^{m\times m}   assumed unit upper triangular
*     V_2 \in \R^{m\times n-m}
*
*     Compute T = V_1'*T
*
      CALL DTRTRM('Left', 'Lower', 'Transpose', 'Non-unit', 'Unit',
     $         M, ONE, Q, LDQ, Q, LDQ)
*
*     Compute Q = -TV. This means that we need to break apart
*     Our computation in two parts
*
*           |---------------|
*     Q =   | -T*V_1 -T*V_2 |
*           |---------------|
*
*     Q_2 = -T*V_2 (TRMM) but only when necessary
*
      IF (N.GT.M) THEN
         CALL DTRMM('Left', 'Lower', 'No Transpose', 'Non-unit',
     $            M, N-M, NEG_ONE, Q, LDQ, Q(1,M+1), LDQ)
      END IF
*
*     Q_1 = -T*V_1 (Lower-Upper Matrix-Matrix multiplication)
*
      CALL DLUMM('Left', 'Non-unit', 'Unit', M, NEG_ONE, Q, LDQ)
*
*     Q = "I" + Q
*
      J = MIN(M,N)
      DO I = 1, J
         Q(I,I) = Q(I,I) + ONE
      END DO
      END SUBROUTINE
