*> \brief \b CLARFT_UT forms the triangular factor T of a block reflector H = I - vtvH
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CLARFT_UT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarft.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarft.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarft.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLARFT_UT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, STOREV
*       INTEGER            K, LDT, LDV, N
*       ..
*       .. Array Arguments ..
*       COMPLEX         T( LDT, * ), TAU( * ), V( LDV, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLARFT_UT forms the triangular factor T of a complex block reflector H
*> of order n, which is defined as a product of k elementary reflectors.
*>
*> If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*>
*> If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*>
*> If STOREV = 'C', the vector which defines the elementary reflector
*> H(i) is stored in the i-th column of the array V, and
*>
*>    H  =  I - V * T * V**T
*>
*> If STOREV = 'R', the vector which defines the elementary reflector
*> H(i) is stored in the i-th row of the array V, and
*>
*>    H  =  I - V**T * T * V
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Specifies the order in which the elementary reflectors are
*>          multiplied to form the block reflector:
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Specifies how the vectors which define the elementary
*>          reflectors are stored (see also Further Details):
*>          = 'C': columnwise
*>          = 'R': rowwise
*> \endverbatim
*>
*> \param[in] JOBT
*> \verbatim
*>          JOBT is CHARACTER*1
*>          Specifies if we are computing a T compatible with larfb
*>          or not
*>          = '1': Compute a compatible T (requires an inverse)
*>          = '2': Do not compute such a T (must use TRSM version of larfb
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the block reflector H. N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the triangular factor T (= the number of
*>          elementary reflectors). K >= 1.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension
*>                               (LDV,K) if STOREV = 'C'
*>                               (LDV,N) if STOREV = 'R'
*>          The matrix V. See further details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i).
*> \endverbatim
*>
*> \param[out] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,K)
*>          The k by k triangular factor T of the block reflector.
*>          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*>          lower triangular. The rest of the array is not used.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
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
*> \ingroup larft
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored.
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*>                   ( v1  1    )                     (     1 v2 v2 v2 )
*>                   ( v1 v2  1 )                     (        1 v3 v3 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*>                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*>                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*>                   (     1 v3 )
*>                   (        1 )
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CLARFT_UT( DIRECT, STOREV, JOBT, N, K, V, LDV, TAU,
     $            T, LDT )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     This subroutine compute the T factor following the algorithm
*     laid out in https://www.cs.utexas.edu/users/flame/pubs/p169-joffrain.pdf
*     Which we adapt to the 5 other cases of computing T we are interested in
*
*
*     .. Scalar Arguments
*
      CHARACTER         DIRECT, STOREV, JOBT
      INTEGER           K, LDT, LDV, N, INFO
*     ..
*     .. Array Arguments ..
*
      COMPLEX           T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*     .. Parameters ..
*
      REAL              ZERO, ONE
      PARAMETER(ZERO = 0.0E+0, ONE = 1.0E+0)
*
*     .. Local Scalars ..
*
      INTEGER           I
      LOGICAL           QR,LQ,QL,RQ,LQT,RQT,DIRF,COLV,TDIRF,TCOLV,INVT
*
*     .. External Subroutines ..
*
      EXTERNAL          CHT3RK, CHERK, CTRTRI
*
*     .. External Functions..
*
      LOGICAL           LSAME
      EXTERNAL          LSAME
*
*     .. Intrinsic Functions..
*
      INTRINSIC         CMPLX,CONJG
*
*     Beginning of executable statements
*
*     This method is only viable for non-singular V matrices with
*     non-zero associated tau values. We know for a fact that V is 
*     always non-singular as V is unit triangular, however tau can
*     be 0. Thus, if we detect this case, we bail to the level-2 BLAS
*     implementation, which is known to work in these instances, which
*     will never bail back to this implementation.
*
*     Note: We are not bailing to the standard larft as that subroutine
*     may call this subroutine as a terminating case.
*     We could also just error out by calling XERBLA if this is desired
*
      DO I = 1, K
         IF( TAU(I).EQ.CMPLX(ZERO) ) THEN
*           TODO: error here if jobt == 2
            CALL CLARFT_LVL2(DIRECT, STOREV, N, K, V, LDV, TAU,
     $            T, LDT)
            RETURN
         END IF
      END DO
*
*     If we reach here, then we guarantee the calls to dtrtri will not fail
*
*     Determine what kind of Q we need to compute
*     We assume that if the user doesn't provide 'F' for DIRECT,
*     then they meant to provide 'B' and if they don't provide
*     'C' for STOREV, then they meant to provide 'R'
*
      DIRF = LSAME(DIRECT,'F')
      TDIRF = LSAME(DIRECT,'T')
      COLV = LSAME(STOREV,'C')
      TCOLV = LSAME(STOREV,'T')
      INVT = LSAME(JOBT, '1')
*
*     QR happens when we have forward direction in column storage
*
      QR = DIRF.AND.COLV
*
*     LQT happens when we have forward direction in row storage and want to
*     compute the conjugate transpose of the T we would normally compute
*
      LQT = DIRF.AND.TCOLV
*
*     LQ happens when we have forward direction in row storage and want to
*     compute the T we would normally compute
*
      LQ = DIRF.AND.(.NOT.LQT)
*
*     QL happens when we have backward direction in column storage
*
      QL = (.NOT.DIRF).AND.COLV
*
*     RQT happens when we have backward direction in row storage and want to
*     compute the conjugate transpose of the T we would normally compute
*
      RQT = TDIRF.AND.(.NOT.COLV)
*
*     RQ happens when we have backward direction in row storage and want to 
*     compute the T that we would normally compute
*
      RQ = (.NOT.RQT).AND.(.NOT.COLV)
*
*     Note that in the following cases, we use ut(A) and lt(A) to
*     denote the upper and lower triangular components of the matrix A
*
      IF(QR) THEN
*
*        Break V apart into 2 components
*
*        V = |-----|
*            | V_1 | k   
*            | V_2 | n-k 
*            |-----|
*              k
*
*        Where V_1 is unit lower triangular and V_2 is a general
*        rectangular matrix
*
*        We compute T as follows:
*
*        T = ut(V**H * V) = ut(V_1**H * V_1 + V_2**H * V_2)
*        diag(T) = tau
*        T = T^{-1}
*
*        Compute T = ut(V_1**H * V_1)
*
         CALL CHT3RK('Lower', 'Upper', 'Conjugate', 'Unit', K, ONE,
     $         V, LDV, ZERO, T, LDT)
*
*        Compute T = ut(V_2**H * V_2 + T)
*
         CALL CHERK('Upper', 'Conjugate', K, N-K, ONE, V(K+1,1),
     $         LDV, ONE, T, LDT)
*
*        Set the diagonals to 1/tau since we first construct T^{-1}
*           Note: we ensured all of these values are non-zero above
*
         DO I = 1, K
            T(I,I) = 1/TAU(I)
         END DO
*
*        Compute T = T^{-1}
*           Note that this cannot fail as we already ensured each tau value was
*           non-zero
*
         IF( INVT ) THEN
            CALL CTRTRI('Upper', 'Non-Unit', K, T, LDT, INFO)
         END IF
      ELSE IF(LQ) THEN
*
*        Break V apart into 2 components
*
*        V = |---------|
*            | V_1 V_2 | k
*            |---------|
*              k   n-k
*
*        Where V_1 is unit upper triangular and V_2 is a general
*        rectangular matrix
*
*        We compute T as follows:
*
*        T = ut(V * V**H) = ut(V_1 * V_1**H + V_2 * V_2**H)
*        diag(T) = tau
*        T = T^{-1}
*
*        Compute T = ut(V_1 * V_1**H)
*
         CALL CHT3RK('Upper', 'Upper', 'No Transpose', 'Unit', K,
     $         ONE, V, LDV, ZERO, T, LDT)
*
*        Compute T = ut(V_2 * V_2**H + T)
*
         CALL CHERK('Upper', 'No Transpose', K, N-K, ONE, V(1,K+1),
     $         LDV, ONE, T, LDT)
*
*        Set the diagonals to 1/tau since we first construct T^{-1}
*           Note: we ensured all of these values are non-zero above
*
         DO I = 1, K
            T(I,I) = 1/TAU(I)
         END DO
*
*        Compute T = T^{-1}
*           Note that this cannot fail as we already ensured each tau value was
*           non-zero
*
         IF( INVT ) THEN
            CALL CTRTRI('Upper', 'Non-Unit', K, T, LDT, INFO)
         END IF
      ELSE IF(LQT) THEN
*
*        Break V apart into 2 components
*
*        V = |---------|
*            | V_1 V_2 | k
*            |---------|
*              k   n-k
*
*        Where V_1 is unit upper triangular and V_2 is a general
*        rectangular matrix
*
*        We compute T as follows:
*
*        T = lt(V * V**H) = lt(V_1 * V_1**H + V_2 * V_2**H)
*        diag(T) = tau
*        T = T^{-1}
*
*        Compute T = lt(V_1 * V_1**H)
*
         CALL CHT3RK('Upper', 'Lower', 'No Transpose', 'Unit', K,
     $         ONE, V, LDV, ZERO, T, LDT)
*
*        Compute T = lt(V_2 * V_2**H + T)
*
         CALL CHERK('Lower', 'No Transpose', K, N-K, ONE, V(1,K+1),
     $         LDV, ONE, T, LDT)
*
*        Set the diagonals to 1/tau since we first construct T^{-1}
*           Note: we ensured all of these values are non-zero above
*
         DO I = 1, K
            T(I,I) = 1/CONJG(TAU(I))
         END DO
*
*        Compute T = T^{-1}
*           Note that this cannot fail as we already ensured each tau value was
*           non-zero
*
         IF( INVT ) THEN
            CALL CTRTRI('Lower', 'Non-Unit', K, T, LDT, INFO)
         END IF
      ELSE IF(QL) THEN
*
*        Break V apart into 2 components
*
*        V = |-----|
*            | V_2 | n-k   
*            | V_1 | k 
*            |-----|
*              k
*
*        Where V_1 is unit upper triangular and V_2 is a general
*        rectangular matrix
*
*        We compute T as follows:
*
*        T = lt(V**H * V) = lt(V_1**H * V_1 + V_2**H * V_2)
*        diag(T) = tau
*        T = T^{-1}
*
*        Compute T = lt(V_1**H * V_1)
*
         CALL CHT3RK('Upper', 'Lower', 'Conjugate', 'Unit', K, ONE,
     $         V(N-K+1,1), LDV, ZERO, T, LDT)
*
*        Compute T = lt(V_2**H * V_2 + T)
*
         CALL CHERK('Lower', 'Conjugate', K, N-K, ONE, V, LDV,
     $         ONE, T, LDT)
*
*        Set the diagonals to 1/tau since we first construct T^{-1}
*           Note: we ensured all of these values are non-zero above
*
         DO I = 1, K
            T(I,I) = 1/TAU(I)
         END DO
*
*        Compute T = T^{-1}
*           Note that this cannot fail as we already ensured each tau value was
*           non-zero
*
         IF( INVT ) THEN
            CALL CTRTRI('Lower', 'Non-Unit', K, T, LDT, INFO)
         END IF
      ELSE IF(RQ) THEN
*
*        Break V apart into 2 components
*
*        V = |---------|
*            | V_2 V_1 | k
*            |---------|
*              n-k k
*
*        Where V_1 is unit lower triangular and V_2 is a general
*        rectangular matrix
*
*        We compute T as follows:
*
*        T = lt(V * V**H) = lt(V_1 * V_1**H + V_2 * V_2**H)
*        diag(T) = tau
*        T = T^{-1}
*
*        Compute T = lt(V_1 * V_1**H)
*
         CALL CHT3RK('Lower', 'Lower', 'No Transpose', 'Unit', K,
     $         ONE, V(1,N-K+1), LDV, ZERO, T, LDT)
*
*        Compute T = lt(V_2 * V_2**H + T)
*
         CALL CHERK('Lower', 'No Transpose', K, N-K, ONE, V, LDV,
     $         ONE, T, LDT)
*
*        Set the diagonals to 1/tau since we first construct T^{-1}
*           Note: we ensured all of these values are non-zero above
*
         DO I = 1, K
            T(I,I) = 1/TAU(I)
         END DO
*
*        Compute T = T^{-1}
*           Note that this cannot fail as we already ensured each tau value was
*           non-zero
*
         IF( INVT ) THEN
            CALL CTRTRI('Lower', 'Non-Unit', K, T, LDT, INFO)
         END IF
      ELSE IF(RQT) THEN
*
*        Break V apart into 2 components
*
*        V = |---------|
*            | V_2 V_1 | k
*            |---------|
*              n-k k
*
*        Where V_1 is unit lower triangular and V_2 is a general
*        rectangular matrix
*
*        We compute T as follows:
*
*        T = ut(V * V**H) = ut(V_1 * V_1**H + V_2 * V_2**H)
*        diag(T) = tau
*        T = T^{-1}
*
*        Compute T = ut(V_1 * V_1**H)
*
         CALL CHT3RK('Lower', 'Upper', 'No Transpose', 'Unit', K,
     $         ONE, V(1,N-K+1), LDV, ZERO, T, LDT)
*
*        Compute T = ut(V_2 * V_2**H + T)
*
         CALL CHERK('Upper', 'No Transpose', K, N-K, ONE, V, LDV,
     $         ONE, T, LDT)
*
*        Set the diagonals to 1/tau since we first construct T^{-1}
*           Note: we ensured all of these values are non-zero above
*
         DO I = 1, K
            T(I,I) = 1/CONJG(TAU(I))
         END DO
*
*        Compute T = T^{-1}
*           Note that this cannot fail as we already ensured each tau value was
*           non-zero
*
         IF( INVT ) THEN
            CALL CTRTRI('Upper', 'Non-Unit', K, T, LDT, INFO)
         END IF
      END IF
      END SUBROUTINE
