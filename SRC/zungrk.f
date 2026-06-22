*> \brief \b ZUNGRK computes the explicit Q factor from DGERQF and ZLARFT
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE ZUNGRK(M, N, Q, LDQ)
*
*        .. Scalar Arguments ..
*        INTEGER           M, N, LDQ
*        ..
*        .. Array Arguments ..
*        COMPLEX*16        Q(LDQ,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZUNGRK generates an m by n complex matrix Q with orthonormal rows,
*> which is defined as the last m rows of the product of m
*> elementary reflectors
*>
*>       Q  =  I - V'*T*V = H(m) . . . H(2) H(1)
*>
*> Where V is an m by n matrix whose columns are householder reflectors
*> as returned by ZGERQF and T is the n by n matrix returned by ZLARFT
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix V, and the order of T.
*>          M >= 0.
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
*>       Q is COMPLEX array, dimension (LDQ,N)
*>       On entry, Q(i,1:n-m-1+i) contains the vector which defines the
*>       elementary reflector H(i), for i=1,...,n as returned by ZGERKF.
*>       In addition, the upper triangular portion of the submatrix given
*>       by Q(1:m,n-m:n) will contain the array T as returned by ZLARFT.
*>       See further details for more information.
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
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*> The storage of the V and T components inside Q is best illustrated by
*> the following example with m = 3, n = 5.
*>
*> Q =   |----------------|
*>       | V1 V1 T1 T1 T1 |
*>       | V2 V2 V2 T2 T2 |
*>       | V3 V3 V3 V3 T3 |
*>       |----------------|
*>
*> \endverbatim
*>
*  =====================================================================

      SUBROUTINE ZUNGRK(M, N, Q, LDQ)
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER           M, N, LDQ
*     ..
*     .. Array Arguments ..
      COMPLEX*16        Q(LDQ,*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16        NEG_ONE, ONE
      PARAMETER(NEG_ONE=(-1.0D+0,0.0D+0), ONE=(1.0D+0,0.0D+0))
*     ..
*     .. Local Scalars ..
      INTEGER           I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL          ZTRMM, ZTRTRM, ZLUMM
*     ..
*     .. Intrinsic Functions..
      INTRINSIC         MIN
*     ..
*     .. Executable Statements ..
*
*     Break Q apart as follows
*
*           |-----|
*     Q =   | V T |
*           |-----|
*
*     Where T is an m-by-m upper triangular matrix, and V is as described
*     in the Further Details section
*
*     In turn, break apart V as follows
*
*           |---------|
*     V =   | V_2 V_1 |
*           |---------|
*
*     Where:
*
*     V_1 \in \R^{m\times m}   assumed unit lower triangular
*     V_2 \in \R^{m\times n-m}
*
*     Compute T = V_1'*T
*
      CALL ZTRTRM('Left', 'Upper', 'Conjugate Transpose',
     $         'Non-Unit', 'Unit', M, ONE, Q(1,N-M+1), LDQ, Q(1,N-M+1),
     &         LDQ)
*
*     Compute Q = -TV. This means that we need to break apart
*     Our computation in two parts
*
*           |---------------|
*     Q =   | -T*V_2 -T*V_1 |
*           |---------------|
*
*     Q_2 = -T*V_2 (TRMM) but only when necessary
*
      IF (N.GT.M) THEN
         CALL ZTRMM('Left', 'Upper', 'No Transpose', 'Non-Unit',
     $            M, N-M, NEG_ONE, Q(1,N-M+1), LDQ, Q, LDQ)
      END IF
*
*     Q_1 = -T*V_1 (Lower-Upper Matrix-Matrix multiplication)
*
      CALL ZLUMM('Right', 'Unit', 'Non-Unit', M, NEG_ONE,
     $         Q(1,N-M+1), LDQ)
*
*     Q = "I" + Q
*
      J = MIN(M,N)
      DO I = 1, J
         Q(I,N-M+I) = Q(I,N-M+I) + ONE
      END DO
      END SUBROUTINE
