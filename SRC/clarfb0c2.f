*> \brief \b CLARFB0C2 applies a block reflector or its conjugate-transpose 
* to a rectangular matrix with a 0 block while constructing the explicit Q
* factor
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*
*  Definition:
*  ===========
*
*     SUBROUTINE CLARFB0C2(C2I, SIDE, TRANS, DIRECT, STOREV, M, N,
*    $                     K, V, LDV, T, LDT, C, LDC)
*        ! Scalar arguments
*        INTEGER           M, N, K, LDV, LDC, LDT
*        CHARACTER         SIDE, TRANS, DIRECT, STOREV 
*        ! True means that we are assuming C2 is the identity matrix
*        !     and thus don't reference whatever is present in C2 
*        !     at the beginning.
*        LOGICAL           C2I
*        ! Array arguments
*        COMPLEX           V(LDV,*), C(LDC,*), T(LDT,*)
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLARFB0C2 applies a real block reflector H or its transpose H**H to a
*> complex m by n matrix C with a 0 block, while computing the explicit Q factor
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] C2I
*> \verbatim
*>          C2I is LOGICAL
*>          = .TRUE.: Assume the nonzero block of C is the identity matrix
*>          = .FALSE.: Use existing data in the nonzero block of C
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply H or H**H from the Left
*>          = 'R': apply H or H**H from the Right
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          = 'N': apply H (No transpose)
*>          = 'C': apply H**H (Conjugate transpose)
*> \endverbatim
*>
*> \param[in] DIRECT
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          Indicates how H is formed from a product of elementary
*>          reflectors
*>          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*>          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*> \endverbatim
*>
*> \param[in] STOREV
*> \verbatim
*>          STOREV is CHARACTER*1
*>          Indicates how the vectors which define the elementary
*>          reflectors are stored:
*>          = 'C': Columnwise
*>          = 'R': Rowwise
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the matrix T (= the number of elementary
*>          reflectors whose product defines the block reflector).
*>          If SIDE = 'L', M >= K >= 0;
*>          if SIDE = 'R', N >= K >= 0.
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension
*>                                (LDV,K) if STOREV = 'C'
*>                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*>                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*>          See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*>          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*>          if STOREV = 'R', LDV >= K.
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,K)
*>          The triangular K-by-K matrix T in the representation of the
*>          block reflector.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension (LDC,N)
*>          On entry, the M-by-N matrix C.
*>          On exit, C is overwritten by H*C or H**H*C or C*H or C*H**H.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
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
*> \ingroup larfb
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The triangular part of V (including its diagonal) is not
*>  referenced.
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
      SUBROUTINE CLARFB0C2(C2I, SIDE, TRANS, DIRECT, STOREV,
     $                     M, N, K, V, LDV, T, LDT, C, LDC)
         ! Scalar arguments
         INTEGER           M, N, K, LDV, LDC, LDT
         CHARACTER         SIDE, TRANS, DIRECT, STOREV 
         ! True means that we are assuming C2 is the identity matrix
         !     and thus don't reference whatever is present in C2 
         !     at the beginning.
         LOGICAL           C2I

         ! Array arguments
         COMPLEX           V(LDV,*), C(LDC,*), T(LDT,*)
         ! Local scalars
         LOGICAL           QR, LQ, QL, DIRF, COLV, SIDEL, SIDER,
     $                     TRANST
         INTEGER           I, J
         ! Intrinsic Functions
         INTRINSIC         CONJG
         ! External functions
         LOGICAL           LSAME
         EXTERNAL          LSAME
         ! External Subroutines
         EXTERNAL          CGEMM, CTRMM
         ! Parameters
         COMPLEX           ONE, ZERO, NEG_ONE
         PARAMETER(ONE=(1.0E+0, 0.0E+0),
     $            ZERO = (0.0E+0, 0.0E+0), 
     $            NEG_ONE = (-1.0E+0, 0.0E+0))

         ! Beginning of executable statements
         ! Convert our character flags to logical values
         DIRF = LSAME(DIRECT,'F')
         COLV = LSAME(STOREV,'C')
         SIDEL = LSAME(SIDE,'L')
         SIDER = LSAME(SIDE,'R')
         TRANST = LSAME(TRANS,'C')

         ! Determine which of the 4 modes are using.
         ! QR is when we store the reflectors column by column and have the
         ! 'first' reflector stored in the first column
         QR = DIRF.AND.COLV

         ! LQ is when we store the reflectors row by row and have the
         ! 'first' reflector stored in the first row
         LQ = DIRF.AND.(.NOT.COLV)

         ! QL is when we store the reflectors column by column and have the
         ! 'first' reflector stored in the last column
         QL = (.NOT.DIRF).AND.COLV

         ! RQ is when we store the reflectors row by row and have the
         ! 'first' reflector stored in the last row
         ! RQ = (.NOT.DIRF).AND.(.NOT.COLV)
         ! Since we have exactly one of these 4 modes, we don't need to actually
         ! store the value of RQ, instead we assume this is the case if we fail
         ! the above 3 checks.

         IF (QR) THEN
            ! We are computing C = HC = (I - VTV')C
            ! Where: V = [ V1 ] and C = [ C1 ]
            !            [ V2 ]         [ C2 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{M-K\times K}
            !     C1=0\in\R^{K\times N}
            !     C2\in\R^{M-K\times N}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = HC = (I - [ V1 ]T [V1', V2'])[ 0  ]
            !               [ V2 ]             [ C2 ]
            !   = [ 0  ] - [ V1 ]T*V2'*C2
            !     [ C2 ]   [ V2 ]
            !
            !   = [ 0  ] - [ V1*T*V2'*C2 ]
            !     [ C2 ]   [ V2*T*V2'*C2 ]
            !
            !   = [      V1*T*V2'*C2 ]
            !     [ C2 - V2*T*V2'*C2 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = V2'*C2
            ! C1 = T*C1
            ! C2 = C2 - V2*C1
            ! C1 = -V1*C1
            !
            ! To achieve the same end result
            !
            ! Check to ensure side and trans are the expected values 
            !
            IF( .NOT.SIDEL ) THEN
               CALL XERBLA('CLARFB0C2', 2)
               RETURN
            ELSE IF(TRANST) THEN
               CALL XERBLA('CLARFB0C2', 3)
               RETURN
            END IF
            !
            ! C1 = V2'*C2
            !
            IF (C2I) THEN
               DO J = 1, N
                  DO I = 1, K
                     C(I,J) = CONJG(V(K+J,I))
                  END DO
               END DO
            ELSE
               CALL CGEMM('Conjugate', 'No Transpose', K, N, M - K,
     $                     ONE, V(K+1,1), LDV, C(K+1,1), LDC, ZERO,
     $                     C, LDC)
            END IF
            !
            ! C1 = T*C1
            !
            CALL CTRMM('Left', 'Upper', 'No Transpose', 'Non-unit',
     $                  K, N, ONE, T, LDT, C, LDC)
            !
            ! C2 = C2 - V2*C1 = -V2*C1 + C2
            !
            IF (C2I) THEN
               CALL CGEMM('No Transpose', 'No Transpose', M-K, N, K,
     $                     NEG_ONE, V(K+1,1), LDV, C, LDC, ZERO,
     $                     C(K+1,1), LDC)
               DO I = 1, N
                  C(K+I,I) = C(K+I,I) + ONE
               END DO
            ELSE
               CALL CGEMM('No Transpose', 'No Transpose', M-K, N, K,
     $                     NEG_ONE, V(K+1,1), LDV, C, LDC, ONE,
     $                     C(K+1,1), LDC)
            END IF
            !
            ! C1 = -V1*C1
            !
            CALL CTRMM('Left', 'Lower', 'No Transpose', 'Unit',
     $                  K, N, NEG_ONE, V, LDV, C, LDC)
         ELSE IF (LQ) THEN
            ! We are computing C = CH' = C(I-V'T'V)
            ! Where: V = [ V1 V2 ] and C = [ C1 C2 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{K\times N-K}
            !     C1=0\in\R^{M\times K}
            !     C2\in\R^{M\times N-K}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = CH' = [ 0, C2 ](I - [ V1' ]T'[ V1, V2 ])
            !                         [ V2' ]
            !
            !   = [ 0, C2 ] - [ 0, C2 ][ V1' ]T'[ V1, V2 ]
            !                          [ V2' ]
            !
            !   = [ 0, C2 ] - C2*V2'*T'[ V1, V2 ]
            !
            !   = [ -C2*V2'*T'*V1, C2 - C2*V2'*T'*V2 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = C2*V2'
            ! C1 = C1*T'
            ! C2 = C2 - C1*V2
            ! C1 = -C1*V1
            !
            ! To achieve the same end result
            !
            ! Check to ensure side and trans are the expected values 
            !
            IF( .NOT.SIDER ) THEN
               CALL XERBLA('CLARFB0C2', 2)
               RETURN
            ELSE IF(.NOT.TRANST) THEN
               CALL XERBLA('CLARFB0C2', 3)
               RETURN
            END IF
            !
            ! C1 = C2*V2'
            !
            IF( C2I ) THEN
               DO J = 1, K
                  DO I = 1, M
                     C(I,J) = CONJG(V(J,K+I))
                  END DO
               END DO
            ELSE
               CALL CGEMM('No Transpose', 'Conjugate', M, K, N-K,
     $               ONE, C(1,K+1), LDC, V(1, K+1), LDV, ZERO, C,
     $               LDC)
            END IF
            !
            ! C1 = C1*T'
            !
            CALL CTRMM('Right', 'Upper', 'Conjugate', 'Non-unit',
     $            M, K, ONE, T, LDT, C, LDC)
            !
            ! C2 = C2 - C1*V2 = -C1*V2 + C2
            !
            IF( C2I ) THEN
               CALL CGEMM('No Transpose', 'No Transpose', M, N-K, K,
     $               NEG_ONE, C, LDC, V(1,K+1), LDV, ZERO, C(1,K+1),
     $               LDC)
               DO I = 1, M
                  C(I,K+I) = C(I,K+I) + ONE
               END DO
            ELSE
               CALL CGEMM('No Transpose', 'No Transpose', M, N-K, K,
     $               NEG_ONE, C, LDC, V(1,K+1), LDV, ONE, C(1,K+1),
     $               LDC)
            END IF
            !
            ! C1 = -C1*V1
            !
            CALL CTRMM('Right', 'Upper', 'No Transpose', 'Unit',
     $            M, K, NEG_ONE, V, LDV, C, LDC)
         ELSE IF (QL) THEN
            ! We are computing C = HC = (I - VTV')C
            ! Where: V = [ V2 ] and C = [ C2 ]
            !            [ V1 ]         [ C1 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{M-K\times K}
            !     C1=0\in\R^{K\times N}
            !     C2\in\R^{M-K\times N}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = HC = (I-[ V2 ]T[ V2' V1' ])[ C2 ]
            !             [ V1 ]             [ 0  ]
            !
            !   = [ C2 ] - [ V2 ]T*V2'*C2
            !     [ 0  ]   [ V1 ]
            !
            !   = [ C2 ] - [ V2*T*V2'*C2 ]
            !     [ 0  ]   [ V1*T*V2'*C2 ]
            !
            !   = [ C2 - V2*T*V2'*C2 ]
            !     [    - V1*T*V2'*C2 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = V2'*C2
            ! C1 = T*C1
            ! C2 = C2 - V2*C1
            ! C1 = -V1*C1
            !
            ! To achieve the same end result
            !
            ! Check to ensure side and trans are the expected values 
            !
            IF( .NOT.SIDEL ) THEN
               CALL XERBLA('CLARFB0C2', 2)
               RETURN
            ELSE IF(TRANST) THEN
               CALL XERBLA('CLARFB0C2', 3)
               RETURN
            END IF
            !
            ! C1 = V2'*C2
            !
            IF( C2I ) THEN
               DO J = 1, N
                  DO I = 1, K
                     C(M-K+I,J) = CONJG(V(J,I))
                  END DO
               END DO
            ELSE
               CALL CGEMM('Conjugate', 'No Transpose', K, N, M-K,
     $            ONE, V, LDV, C, LDC, ZERO, C(M-K+1, 1), LDC)
            END IF
            !
            ! C1 = T*C1
            !
            CALL CTRMM('Left', 'Lower', 'No Transpose', 'Non-unit',
     $         K, N, ONE, T, LDT, C(M-K+1,1), LDC)
            !
            ! C2 = C2 - V2*C1 = -V2*C1 + C2
            !
            IF( C2I ) THEN
               CALL CGEMM('No Transpose', 'No Transpose', M-K, N, K,
     $            NEG_ONE, V, LDV, C(M-K+1,1), LDC, ZERO, C, LDC)
               DO I = 1, N
                  C(I,I) = C(I,I) + ONE
               END DO
            ELSE
               CALL CGEMM('No Transpose', 'No Transpose', M-K, N, K,
     $            NEG_ONE, V, LDV, C(M-K+1,1), LDC, ONE, C, LDC)
            END IF
            !
            ! C1 = -V1*C1
            !
            CALL CTRMM('Left', 'Upper', 'No Transpose', 'Unit',
     $         K, N, NEG_ONE, V(M-K+1,1), LDV, C(M-K+1,1), LDC)
         ELSE ! IF (RQ) THEN
            ! We are computing C = CH' = C(I-V'T'V)
            ! Where: V = [ V2 V1] and C = [ C2 C1 ]
            ! with the following dimensions:
            !     V1\in\R^{K\times K}
            !     V2\in\R^{K\times N-K}
            !     C1=0\in\R^{M\times K}
            !     C2\in\R^{M\times N-K}
            ! Since we are assuming that C1 is a zero matrix and it will be
            ! overwritten on exit, we can use this spot as a temporary workspace
            ! without having to allocate anything extra.
            ! This lets us simplify our above equation to get
            !
            ! C = CH' = [ C2, 0 ] (I - [ V2' ]T'[ V2, V1 ]
            !                          [ V1' ]
            !
            !   = [ C2, 0 ] - [ C2, 0 ] [ V2' ]T'[ V2, V1 ]
            !                           [ V1' ]
            !
            !   = [ C2, 0 ] - C2*V2'*T'[ V2, V1 ]
            !
            !   = [ C2, 0 ] - [ C2*V2'*T'*V2, C2*V2'*T'*V1 ]
            !
            !   = [ C2 - C2*V2'*T'*V2, -C2*V2'*T'*V1 ]
            !
            ! So, we can order our computations as follows:
            !
            ! C1 = C2*V2'
            ! C1 = C1*T'
            ! C2 = C2 - C1*V2
            ! C1 = -C1*V1
            !
            !
            ! To achieve the same end result
            !
            ! Check to ensure side and trans are the expected values 
            !
            IF( .NOT.SIDER ) THEN
               CALL XERBLA('CLARFB0C2', 2)
               RETURN
            ELSE IF(.NOT.TRANST) THEN
               CALL XERBLA('CLARFB0C2', 3)
               RETURN
            END IF
            !
            ! C1 = C2*V2'
            !
            IF( C2I ) THEN
               DO J = 1, K
                  DO I = 1, M
                     C(I,N-K+J) = CONJG(V(J,I))
                  END DO
               END DO
            ELSE
               CALL CGEMM('No Transpose', 'Conjugate', M, K, N-K,
     $            ONE, C, LDC, V, LDV, ZERO, C(1, N-K+1), LDC)
            END IF
            !
            ! C1 = C1*T'
            !
            CALL CTRMM('Right', 'Lower', 'Conjugate', 'Non-unit',
     $         M, K, ONE, T, LDT, C(1, N-K+1), LDC)
            !
            ! C2 = C2 - C1*V2 = -C1*V2 + C2
            !
            IF( C2I ) THEN
               CALL CGEMM('No Transpose', 'No Transpose', M, N-K, K,
     $            NEG_ONE, C(1, N-K+1), LDC, V, LDV, ZERO, C, LDC)
               DO I = 1, M
                  C(I,I) = C(I,I) + ONE
               END DO
            ELSE
               CALL CGEMM('No Transpose', 'No Transpose', M, N-K, K,
     $            NEG_ONE, C(1, N-K+1), LDC, V, LDV, ONE, C, LDC)
            END IF
            !
            ! C1 = -C1*V1
            !
            CALL CTRMM('Right', 'Lower', 'No Transpose', 'Unit',
     $         M, K, NEG_ONE, V(1, N-K+1), LDV, C(1,N-K+1), LDC)
         END IF
      END SUBROUTINE
