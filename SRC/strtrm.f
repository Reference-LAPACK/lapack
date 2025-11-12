*> \brief \b STRTRM computes an in place triangular-triangular matrix product
*
*  =========== DOCUMENTATION ===========
*
*  Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     RECURSIVE SUBROUTINE STRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV,
*    $                        N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
*        INTEGER           N, LDT, LDV
*        CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
*        REAL              ALPHA
*        ..
*        .. Array Arguments ..
*        REAL              T(LDT,*), V(LDV,*)
*        ..
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> STRTRM performs one  of the matrix-matrix operations
*>
*>       T = \alpha op(V) * T
*>                      or
*>       T = \alpha T * op(V)
*> where \alpha is a scalar, T and V are unit, or non-unit, upper or
*> lower triangular matrix, and op(V) is one of
*>
*>       op(V) = V      or       op(V) = V**T
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>           On entry, SIDE specifies whether op(V) multiplies T from
*>           the left or right as follows:
*>
*>             SIDE = 'L' or 'l'    T = \alpha op(V) * T
*>
*>             SIDE = 'R' or 'r'    T = \alpha T * op(V)
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix T is an upper or
*>           lower triangular matrix as follows:
*>             UPLO = 'U' or 'u'    T is upper triangular
*>
*>             UPLO = 'L' or 'l'    T is lower triangular
*> \Endverbatim
*>
*> \param[in] TRANSV
*> \verbatim
*>          TRANSV is CHARACTER*1
*>           On entry, TRANSV specifies the form of op(V) to be used in
*>           the matrix multiplication as follows:
*>             TRANSV = 'N' or 'n'    op(V) = V
*>
*>             TRANSV = 'T' or 't'    op(V) = V**T
*>
*>             TRANSV = 'C' or 'c'    op(V) = V**T
*> \endverbatim
*>
*> \param[in] DIAGT
*> \verbatim
*>          DIAGT is CHARACTER*1
*>           On entry, DIAGT specifies whether or not T is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      T is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      T is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] DIAGV
*> \verbatim
*>          DIAGV is CHARACTER*1
*>           On entry, DIAGV specifies whether or not V is unit triangular
*>           as follows:
*>
*>              DIAG = 'U' or 'u'      V is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'      V is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of rows and columns of T.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is REAL            .
*>           On entry, ALPHA specifies the scalar alpha. When alpha is
*>           zero then T and V are not referenced, and T and V need not
*>           be set before entry.
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is REAL array, dimension ( LDT, N )
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array T must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           T is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array T must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           T is not referenced.
*>           Note that when  DIAGT = 'U' or 'u',  the diagonal elements of
*>           T  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>           On entry, LDT specifies the first dimension of T as declared
*>           in the calling (sub) program. LDT must be at least max( 1, n ).
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is REAL array, dimension ( LDV, N )
*>           Before entry with UPLO = 'U' or 'u', the leading k-by-k
*>           upper triangular part of the array op(V) must contain the upper
*>           triangular matrix and the strictly lower triangular part of
*>           V is not referenced.
*>           Before entry  with  UPLO = 'L' or 'l', the leading k-by-k
*>           lower triangular part of the array op(V) must contain the lower
*>           triangular matrix and the strictly upper triangular part of
*>           V is not referenced.
*>           Note that when  DIAGV = 'U' or 'u',  the diagonal elements of
*>           V  are not referenced either,  but are assumed to be  unity.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>           On entry, LDV specifies the first dimension of T as declared
*>           in the calling (sub) program. LDV must be at least max( 1, n ).
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
      RECURSIVE SUBROUTINE STRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV,
     $                        N, ALPHA, T, LDT, V, LDV)
*
*        .. Scalar Arguments ..
         INTEGER           N, LDT, LDV
         CHARACTER         SIDE, UPLO, TRANSV, DIAGT, DIAGV
         REAL              ALPHA
*        ..
*        .. Array Arguments ..
         REAL              T(LDT,*), V(LDV,*)
*        ..
*
*  =====================================================================
*
*        .. External Functions ..
         LOGICAL           LSAME
         EXTERNAL          LSAME
*        ..
*        .. External Subroutines ..
         EXTERNAL          STRMM, STRMMOOP, SLASET
*        ..
*        .. Local Scalars ..
         INTEGER           K
         LOGICAL           TLEFT, TUPPER, VTRANS, VUNIT, TUNIT
*        ..
*        .. Local Parameters ..
         REAL             ONE, ZERO
         PARAMETER(ONE=1.0E+0, ZERO=0.0E+0)
*        ..
*
*        Beginning of Executable Statements
*
*
*        Early Termination Criteria
*
         IF (ALPHA.EQ.ZERO) THEN
*
*           If ALPHA is 0, then we are just setting T to be the 0 matrix
*
            CALL SLASET(UPLO, N, N, ZERO, ZERO, T, LDT)
            RETURN
         END IF
         TUNIT = LSAME(DIAGT, 'U')
         VUNIT = LSAME(DIAGV, 'U')
*
*        Terminating Case
*
         IF (N.EQ.1) THEN
            IF (VUNIT.AND.TUNIT) THEN
               T(1,1) = ALPHA
            ELSE IF (VUNIT) THEN
               T(1,1) = ALPHA*T(1,1)
            ELSE IF (TUNIT) THEN
               T(1,1) = ALPHA*V(1,1)
            ELSE
               T(1,1) = ALPHA*T(1,1)*V(1,1)
            END IF
            RETURN
         ELSE IF(N.LE.0) THEN
            RETURN
         END IF
*
*        Recursive case
*
         TUPPER = LSAME(UPLO,'U')
         TLEFT  = LSAME(SIDE,'R')
         VTRANS = LSAME(TRANSV,'T').OR.LSAME(TRANSV,'C')

         K = N / 2
         IF(TUPPER) THEN
*
*           T is upper triangular
*
            IF(TLEFT) THEN
*
*              Compute T = T*op(V)
*
               IF(VTRANS) THEN
*
*                 We are computing T = T*V**T, which we break down as follows
*                 |--------------|           |--------------|   |--------------------|
*                 |T_{11}  T_{12}|           |T_{11}  T_{12}|   |V_{11}**T  V_{21}**T|
*                 |0       T_{22}| = \alpha  |0       T_{22}| * |0          V_{22}**T|
*                 |--------------|           |--------------|   |--------------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 V_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha T_{11}*V_{11}**T
*                 T_{12} = \alpha T_{11}*V_{21}**T + \alpha T_{12}*V_{22}**T
*                 T_{22} = \alpha T_{22}*V_{22}**T
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{12} as follows
*
*                 T_{12} = \alpha T_{12}*V_{22}**T          (STRMM)
*                 T_{12} = \alpha T_{11}*V_{21}**T + T_{12} (STRMMOOP)
*
*                 T_{12} = \alpha T_{12}*V_{22}**T
*
                  CALL STRMM('Right', 'Lower', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V(K+1, K+1), LDV, T(1, K+1), LDT)
*
*                 T_{12} = \alpha T_{11}*V_{21}**T + T_{12}
*
                  CALL STRMMOOP('Left', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, K, N-K, ALPHA, T, LDT,
     $                     V(K+1, 1), LDV, ONE, T(1, K+1), LDT)
               ELSE
*
*                 We are computing T = T*V, which we break down as follows
*                 |--------------|           |--------------|   |-------------|
*                 |T_{11}  T_{12}|           |T_{11}  T_{12}|   |V_{11} V_{12}|
*                 |0       T_{22}| = \alpha  |0       T_{22}| * |0      V_{22}|
*                 |--------------|           |--------------|   |-------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha T_{11}*V_{11}
*                 T_{12} = \alpha T_{11}*V_{12} + \alpha T_{12}*V_{22}
*                 T_{22} = \alpha T_{22}*V_{22}
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{12} as follows
*
*                 T_{12} = \alpha T_{12}*V_{22}          (STRMM)
*                 T_{12} = \alpha T_{11}*V_{12} + T_{12} (STRMMOOP)
*
*                 T_{12} = \alpha T_{12}*V_{22}
*
                  CALL STRMM('Right', 'Upper', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V(K+1, K+1), LDV, T(1, K+1), LDT)
*
*                 T_{12} = \alpha T_{11}*V_{21}**T + T_{12}
*
                  CALL STRMMOOP('Left', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, K, N-K, ALPHA, T, LDT,
     $                     V(1, K+1), LDV, ONE, T(1, K+1), LDT)
               END IF
            ELSE
*
*              Compute T = op(V)*T
*
               IF(VTRANS) THEN
*
*                 We are computing T = V**T*T, which we break down as follows
*                 |--------------|           |--------------------|   |--------------|
*                 |T_{11}  T_{12}|           |V_{11}**T  V_{21}**T|   |T_{11}  T_{12}|
*                 |0       T_{22}| = \alpha  |0          V_{22}**T| * |0       T_{22}|
*                 |--------------|           |--------------------|   |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 V_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha V_{11}**T*T_{11}
*                 T_{12} = \alpha V_{11}**T*T_{12} + \alpha V_{21}**T*T_{22}
*                 T_{22} = \alpha V_{22}**T*T_{22}
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{12} as follows
*
*                 T_{12} = \alpha V_{11}**T*T_{12}          (STRMM)
*                 T_{12} = \alpha V_{21}**T*T_{22} + T_{12} (STRMMOOP)
*
*                 T_{12} = \alpha V_{11}**T*T_{12}
*
                  CALL STRMM('Left', 'Lower', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V, LDV, T(1, K+1), LDT)
*
*                 T_{12} = \alpha V_{21}**T*T_{22} + T_{12}
*
                  CALL STRMMOOP('Right', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, K, N-K, ALPHA, T(K+1, K+1),
     $                     LDT, V(K+1, 1), LDV, ONE, T(1, K+1), LDT)
               ELSE
*
*                 We are computing T = V*T, which we break down as follows
*                 |--------------|           |--------------|   |--------------|
*                 |T_{11}  T_{12}|           |V_{11}  V_{12}|   |T_{11}  T_{12}|
*                 |0       T_{22}| = \alpha  |0       V_{22}| * |0       T_{22}|
*                 |--------------|           |--------------|   |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}    T_{12}\in\R^{k\times n-k}
*                                            T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha V_{11}*T_{11}
*                 T_{12} = \alpha V_{11}*T_{12} + \alpha V_{12}*T_{22}
*                 T_{22} = \alpha V_{22}*T_{22}
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{12} as follows
*
*                 T_{12} = \alpha V_{11}*T_{12}          (STRMM)
*                 T_{12} = \alpha V_{12}*T_{22} + T_{12} (STRMMOOP)
*
*                 T_{12} = \alpha V_{11}*T_{12}
*
                  CALL STRMM('Left', 'Upper', TRANSV, DIAGV, K,
     $                     N-K, ALPHA, V, LDV, T(1, K+1), LDT)
*
*                 T_{12} = \alpha V_{12}*T_{22} + T_{12} (STRMMOOP)
*
                  CALL STRMMOOP('Right', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, K, N-K, ALPHA, T(K+1, K+1),
     $                     LDT, V(1, K+1), LDV, ONE, T(1, K+1), LDT)
               END IF
            END IF
         ELSE
*
*           T is lower triangular
*
            IF(TLEFT) THEN
*
*              Compute T = T*op(V)
*
               IF(VTRANS) THEN
*
*                 We are computing T = T*V**T, which we break down as follows
*                 |--------------|           |--------------|   |--------------------|
*                 |T_{11}  0     |           |T_{11}  0     |   |V_{11}**T  0        |
*                 |T_{21}  T_{22}| = \alpha  |T_{21}  T_{22}| * |V_{12}**T  V_{22}**T|
*                 |--------------|           |--------------|   |--------------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha T_{11}*V_{11}**T
*                 T_{21} = \alpha T_{21}*V_{11}**T + \alpha T_{22}*V_{12}**T
*                 T_{22} = \alpha T_{22}*V_{22}**T
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{21} as follows
*
*                 T_{21} = \alpha T_{21}*V_{11}**T          (STRMM)
*                 T_{21} = \alpha T_{22}*V_{12}**T + T_{21} (STRMMOOP)
*
*                 T_{21} = \alpha T_{21}*V_{11}**T
*
                  CALL STRMM('Right', 'Upper', TRANSV, DIAGV, N-K,
     $                     K, ALPHA, V, LDV, T(K+1, 1), LDT)
*
*                 T_{21} = \alpha T_{22}*V_{12}**T + T_{21}
*
                  CALL STRMMOOP('Left', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, N-K, K, ALPHA, T(K+1, K+1),
     $                     LDT, V(1, K+1), LDV, ONE, T(K+1, 1), LDT)
               ELSE
*
*                 We are computing T = T*V, which we break down as follows
*                 |--------------|           |--------------|   |-------------|
*                 |T_{11}  0     |           |T_{11}  0     |   |V_{11} 0     |
*                 |T_{21}  T_{22}| = \alpha  |T_{21}  T_{22}| * |V_{21} V_{22}|
*                 |--------------|           |--------------|   |-------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 V_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha T_{11}*V_{11}
*                 T_{21} = \alpha T_{21}*V_{11} + \alpha T_{22}*V_{21}
*                 T_{22} = \alpha T_{22}*V_{22}
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{21} as follows
*
*                 T_{21} = \alpha T_{21}*V_{11}          (STRMM)
*                 T_{21} = \alpha T_{22}*V_{21} + T_{21} (STRMMOOP)
*
*                 T_{21} = \alpha T_{21}*V_{11}
*
                  CALL STRMM('Right', 'Lower', TRANSV, DIAGV, N-K,
     $                     K, ALPHA, V, LDV, T(K+1, 1), LDT)
*
*                 T_{21} = \alpha T_{22}*V_{12} + T_{21}
*
                  CALL STRMMOOP('Left', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, N-K, K, ALPHA, T(K+1, K+1),
     $                     LDT, V(K+1, 1), LDV, ONE, T(K+1, 1), LDT)
               END IF
            ELSE
*
*              Compute T = op(V)*T
*
               IF(VTRANS) THEN
*
*                 We are computing T = V**T*T, which we break down as follows
*                 |--------------|           |--------------------|   |--------------|
*                 |T_{11}  0     |           |V_{11}**T  0        |   |T_{11}  0     |
*                 |T_{21}  T_{22}| = \alpha  |V_{12}**T  V_{22}**T| * |T_{21}  T_{22}|
*                 |--------------|           |--------------------|   |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}    V_{12}\in\R^{k\times n-k}
*                                            V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha V_{11}**T*T_{11}
*                 T_{21} = \alpha V_{12}**T*T_{11} + \alpha V_{22}**T*T_{21}
*                 T_{22} = \alpha V_{22}**T*T_{22}
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{21} as follows
*
*                 T_{21} = \alpha V_{22}**T*T_{21}          (STRMM)
*                 T_{21} = \alpha V_{12}**T*T_{11} + T_{21} (STRMMOOP)
*
*                 T_{21} = \alpha V_{22}**T*T_{21}
*
                  CALL STRMM('Left', 'Upper', TRANSV, DIAGV, N-K, K,
     $                     ALPHA, V(K+1, K+1), LDV, T(K+1, 1), LDT)
*
*                 T_{21} = \alpha V_{12}**T*T_{11} + T_{21}
*
                  CALL STRMMOOP('Right', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, N-K, K, ALPHA, T, LDT,
     $                     V(1, K+1), LDV, ONE, T(K+1, 1), LDT)
               ELSE
*
*                 We are computing T = V*T, which we break down as follows
*                 |--------------|           |-------------|   |--------------|
*                 |T_{11}  0     |           |V_{11} 0     |   |T_{11}  0     |
*                 |T_{21}  T_{22}| = \alpha  |V_{21} V_{22}| * |T_{21}  T_{22}|
*                 |--------------|           |-------------|   |--------------|
*
*                 Where
*                 T_{11}\in\R^{k\times k}
*                 T_{21}\in\R^{n-k\times k}  T_{22}\in\R^{n-k\times n-k}
*
*                 V_{11}\in\R^{k\times k}
*                 V_{21}\in\R^{n-k\times k}  V_{22}\in\R^{n-k\times n-k}
*
*                 Which means that we get
*
*                 T_{11} = \alpha V_{11}*T_{11}
*                 T_{21} = \alpha V_{21}*T_{11} + \alpha V_{22}*T_{21}
*                 T_{22} = \alpha V_{22}*T_{22}
*
*                 Computing T_{11} and T_{22} are just recursive calls to this
*                 routine, but we can break down computing T_{12} as follows
*
*                 T_{21} = \alpha V_{22}*T_{21}          (STRMM)
*                 T_{21} = \alpha V_{12}*T_{11} + T_{21} (STRMMOOP)
*
*                 T_{21} = \alpha V_{22}*T_{12}
*
                  CALL STRMM('Left', 'Lower', TRANSV, DIAGV, N-K, K,
     $                     ALPHA, V(K+1, K+1), LDV, T(K+1, 1), LDT)
*
*                 T_{21} = \alpha V_{12}*T_{11} + T_{21}
*
                  CALL STRMMOOP('Right', UPLO, 'No Transpose',
     $                     TRANSV, DIAGT, N-K, K, ALPHA, T, LDT,
     $                     V(K+1, 1), LDV, ONE, T(K+1, 1), LDT)
               END IF
            END IF
         END IF
*
*        Since in all the above cases, we compute T_{11} and T_{22}
*        the same, we pass in our flags and call this routine recursively
*
*        Compute T_{11} recursively
*
         CALL STRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV, K, ALPHA,
     $         T, LDT, V, LDV)
*
*        Compute T_{22} recursively
*
         CALL STRTRM(SIDE, UPLO, TRANSV, DIAGT, DIAGV, N-K, ALPHA,
     $         T(K+1, K+1), LDT, V(K+1, K+1), LDV)

      END SUBROUTINE
