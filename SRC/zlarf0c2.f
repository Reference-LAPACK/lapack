*> \brief \b ZLARF0C2 applies an elementary reflector to a rectangular matrix
*> with a 0 row/column while constructing the explicit Q factor.
*
*  =========== DOCUMENTATION ===========
*
* online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  definition:
*  ===========
*
*     SUBROUTINE ZLARF0C2( C2JOB, SIDE, UNITPOS, M, N, TAU,
*                          V, INCV, C, LDC )
*
*     .. scalar arguments ..
*     INTEGER           M, N, LDC
*     CHARACTER         C2JOB, SIDE, UNITPOS
*     COMPLEX  TAU
*     ..
*     .. array arguments ..
*     COMPLEX  C(LDC,*), V( * )
*     ..
*
*
*> \par purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARF0C2 applies a dble elementary reflector H to a dble M-by-N
*> matrix C with a 0 row/column , from either the left or the right.
*> H is represented in the
*> form
*>
*>       H = I - \tau*v*v**T
*>
*> where tau is a dble scalar and v is a dble vector.
*>
*> if tau = 0, then H is taken to be the unit matrix.
*>
*> \endverbatim
*
*  arguments:
*  ==========
*
*> \param[in] C2JOB
*> \verbatim
*>          C2JOB is CHARACTER*1
*>          = 'I': assume C2 is the identity matrix
*>          otherwise: treat C2 as a general matrix that we reference
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply H from the left
*>          = 'R': apply H from the right
*> \endverbatim
*>
*> \param[in] UNITPOS
*> \verbatim
*>          DIRECT is CHARACTER*1
*>          indicates where in v the assumed 1 is located (or if there is one)
*>          = 'F': V has the assumed 1 in the 'first' element of v.
*>          = 'L': V has the assumed 1 in the 'last' element of v.
*>                see further details for more information
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          the number of rows of the matrix C. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          the number of columns of the matrix C. N >= 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX
*>          TAU must contain the value TAU returned from ZLARFG
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          the vector v in the representation of H. V is not used if
*>          TAU = 0. see further details.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          the increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX array, dimension (LDC,N)
*>          on entry, the m by n matrix C.
*>          on exit, C is overwritten by H*C or C*H.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          the leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*
*  authors:
*  ========
*
*> \author univ. of tennessee
*> \author univ. of california berkeley
*> \author univ. of colorado denver
*> \author NAG ltd.
*
*> \par further details:
*  =====================
*>
*> \verbatim
*>
*>  the storage of the vector is best illustrated with
*>  the following example. we use a v of length 5, (this is M if
*>  SIDE = 'L' and N if SIDE = 'R')
*>
*>  UNITPOS = 'F':         UNITPOS = 'L' :
*>
*>          v = (  1 )     v = ( v1 )
*>              ( v2 )         ( v2 )
*>              ( v3 )         ( v3 )
*>              ( v4 )         ( v4 )
*>              ( v5 )         (  1 )
*>
*>  also see zlarf1f and zlarf1l along with their usage in org{2r,2l,r2,l2}
*>  for a similar routine but without the assumed 0 block in C.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZLARF0C2( C2JOB, SIDE, UNITPOS, M, N, TAU,
     $                     V, INCV, C, LDC )
*
*     .. scalar arguments ..
      INTEGER           M, N, LDC, INCV
      CHARACTER         C2JOB, SIDE, UNITPOS
      COMPLEX*16 TAU
*     ..
*     .. array arguments ..
      COMPLEX*16 C(LDC,*), V( * )
*     ..
*
*  =====================================================================
*
*     .. parameters ..
      COMPLEX*16 ONE, ZERO
      PARAMETER         ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. local scalars ..
      LOGICAL           SIDEL, SIDER, C2I, FIRST,
     $                  QR, LQ, QL, RQ
      INTEGER           I, J, LENV
*     ..
*     .. external functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
*     ..
*     .. external subroutines ..
      EXTERNAL          ZGEMM, ZGERU, ZLASET, XERBLA
*     ..
*     .. executable statements ..
*
*     quick return if possible. when \tau is 0, H is defined to be I,
*     so the application doesn't change C at all. however we still flow
*     through as we will need to set the appropriate row of C to the 0 row
*     in addition, our blas kernels should do nothing if the input scalar
*     is 0.
*
      IF( M.LE.0.OR.N.LE.0 ) THEN
         RETURN
      END IF
      FIRST = LSAME(UNITPOS, 'F')
      C2I = LSAME(C2JOB, 'I')
      SIDEL = LSAME(SIDE, 'L')
      SIDER = LSAME(SIDE, 'R')
      IF( .NOT.(SIDEL.OR.SIDER) ) THEN
         CALL XERBLA('ZLARF0C2', 2)
         RETURN
      END IF
*
*
*     determine the length of V for use during an early exit check
*
      LENV = M
      IF( SIDER ) THEN
         LENV = N
      END IF
*
*     if V has length 1, then we just need to set C to be
*     0 and return. this is because we assume that C has 0s
*     in the first/last row/column (the same part that would be touched
*     by the assumed 1 in V), so we should do nothing. however,
*     C may not enter with these 0s set and we need to guarantee
*     these 0s are present on exit. this may need to be changed if
*     we change the above assumption on C's 0 block. note that
*     if lenv is 1, then at least one of m and n is 1, we just don't
*     care in the call to laset
*
      IF( LENV.EQ.1 ) THEN
         CALL ZLASET( 'All', M, N, ZERO, ZERO, C, LDC)
         RETURN
      END IF

*
*     our flags inform what parent routine is calling this. in order
*     to make it's connection to my later blocked version, I will
*     go ahead and do this.
*
      QR = FIRST.AND.SIDEL
      QL = (.NOT.FIRST).AND.SIDEL
      LQ = FIRST.AND.(.NOT.SIDEL)
      RQ = (.NOT.FIRST).AND.(.NOT.SIDEL)
      IF( QR ) THEN
*
*        we are computing C = HC = (I - vtv')C
*        where: v = [ 1 ], C = [ C1 ], and T=tau is a scalar
*                   [ v ]      [ C2 ]
*        with the following dimensions:
*            v\in\R^{m-1}
*            C1=0\in\R^{1\times n}
*            C2  \in\R^{m-1\times n}
*        since we are assuming that C1 is a zero row and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        this lets us simplify our above equation to get
*
*        C = HC = (I - [  1 ]T [ 1 v'])[ 0  ]
*                      [  v ]          [ C2 ]
*          = [ 0  ] - [  1 ]T* v'*C2
*            [ C2 ]   [  v ]
*
*          = [ 0  ] - [ T* v'*C2 ]
*            [ C2 ]   [  v*T* v'*C2 ]
*
*          = [        -T* v'*C2 ]
*            [ C2 -  v*T* v'*C2 ]
*
*        so, we can order our computations as follows:
*
*        C1 = -tau * v'*C2
*        C2 = C2 + v*C1
*
*        if we also add the constraint that C2 starts as the
*        first n columns of the m-1 identity, then this simplifies into
*
*        C1 = -tau*v(1:n)'
*        C2 = I + v*C1
*
*        first check if T = tau is 0, if this is the case, all we
*        need to do is set the first row of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO J = 1, N
               C(1,J) = ZERO
            END DO
*
*           if we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF ( C2I ) THEN
               CALL ZLASET('All', M-1, N, ZERO, ONE,
     $            C(2,1), LDC)
            END IF
         ELSE
            IF( C2I ) THEN
*
*              C1 = -tau*v(1:n)'
*
               DO J = 1, N
                  C(1,J) = -TAU*CONJG(V(1 + (J-1)*INCV))
               END DO
*
*              C2 = I + v*C1
*              no routines that perform this operation exist, so
*              we compute this columnwise
*
               DO J = 1, N
                  DO I = 1, M-1
                     IF (I.EQ.J) THEN
                        C(1+J,J) = ONE + V(1 + (J-1)*INCV)*C(1,J)
                     ELSE
                        C(1+I,J) = V(1 + (I-1)*INCV)*C(1,J)
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau * v'*C2
               CALL ZGEMM('Conjugate Transpose', 'No Transpose',
     $            1, N, M-1, -TAU, V, M, C(2,1), LDC, ZERO, C, LDC)
*
*              C2 = C2 + v*C1
*
               CALL ZGERU(M-1, N, ONE, V, INCV, C, LDC,
     $            C(2,1), LDC)
            END IF
         END IF
      ELSE IF( LQ ) THEN
*        the assumed 1 in V is first in storage: ie our reflector
*        is stored as [ 1 v ]
*
*        we are computing C = CH = C(I - V'TV)
*        where: V = [ 1 v ], C = [ C1 C2 ], and T=tau is a scalar
*
*        with the following dimensions:
*            v   \in\R^{1\times n-1}
*            C1=0\in\R^{m\times   1}
*            C2  \in\R^{m\times n-1}
*        since we are assuming that C1 is a zero column, and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        this lets us simplify our above equation to get
*
*        C = CH = [ 0 C2 ] (I - [ 1  ] T [ 1 v ] )
*                               [ v' ]
*
*               = [ 0 C2 ] - [ 0 C2 ] [ 1  ] T [ 1 v ]
*                                     [ v' ]
*
*               = [ 0 C2 ] - c2v' T [ 1 v ]
*
*               = [ -C2*v'*T, C2 - C2*v'*T*v ]
*
*        so, we can order our computations as follows:
*
*        C1 = -tau * C2*v'
*        C2 = C2 + C1*v
*
*        if we also add the constraint that C2 starts as the
*        first m rows of the n-1 identity, then this simplifies into
*
*        C1 = -tau*v(1:m)'
*        C2 = I + C1*v
*
*        first check if T = tau is 0, if this is the case, all we
*        need to do is set the first column of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO I = 1, M
               C(I,1) = ZERO
            END DO
*
*           if we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF( C2I ) THEN
               CALL ZLASET('All', M, N-1, ZERO, ONE,
     $            C(1,2), LDC)
            END IF
         ELSE
            IF( C2I ) THEN
*
*              C1 = -tau*v(1:m)'
*
               DO I = 1, M
                  C(I,1) = -TAU*CONJG(V(1 + (I-1)*INCV))
               END DO
*
*              C2 = I + C1*v
*
               DO J = 1, N-1
                  DO I = 1, M
                     IF( I.EQ.J ) THEN
                        C(I,J+1) = ONE +
     $                     C(I,1) * V(1 + (J-1) * INCV )
                     ELSE
                        C(I,J+1) = C(I,1) * V(1 + (J-1) * INCV )
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau * C2*v'
*
               CALL ZGEMM('No Transpose', 'Conjugate Transpose',
     $            M, 1, N-1, -TAU, C(1,2), LDC, V, INCV, ZERO,
     $            C, LDC)
*
*              C2 = C2 + C1*v
*
               CALL ZGERU(M, N-1, ONE, C, 1, V, INCV,
     $            C(1,2), LDC)
            END IF
         END IF
      ELSE IF( QL ) THEN
*        the assumed 1 in V is last in storage: ie our refector
*        is stored as [ v ]
*                     [ 1 ]
*        we are computing C = HC = (I - VTV')C
*        where: V = [ v ] and C = [ C2 ] and T=tau is a scalar
*                   [ 1 ]         [ C1 ]
*        with the following dimensions:
*            v\in\R^{m-1\times 1}
*            C1=0\in\R^{1\times n}
*            C2  \in\R^{m-1\times n}
*        since we are assuming that C1 is a zero row and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        this lets us simplify our above equation to get
*
*        C = HC = [ C2 ] - [ v ]T[ v' 1 ][ C2 ]
*                 [ 0  ]   [ 1 ]         [ 0  ]
*
*        = [ C2 ] - [ v ]T*v'C2
*          [ 0  ]   [ 1 ]
*
*        = [ C2 - v*T*v'*C2 ]
*          [    -   T*v'*C2 ]
*
*        so, we can order our computations as follows:
*
*        C1 = -tau*v'*C2
*        C2 = C2 + v*C1
*
*        if we also add the constraint that C2 starts as the
*        last n columns of the m-1 identity, then this simplifies into
*
*        C1 = -tau*v(m-n:m-1)'
*        C2 = I + v*C1
*
*        first check if T = tau is 0, if this is the case, all we
*        need to do is set the last row of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO J = 1, N
               C(M, J) = ZERO
            END DO
*
*           if we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF( C2I ) THEN
               CALL ZLASET('All', M-1, N, ZERO, ZERO, C, LDC)
               DO I = 1, N
                  C(M-N-1+I,I) = ONE
               END DO
            END IF
         ELSE
            IF( C2I ) THEN
*
*              C1 = -tau*v(m-n:m-1)'
*
               DO J = 1, N
                  C(M,J) = -TAU*CONJG(V(1 + (M-N-1 + J-1)*INCV))
               END DO
*
*              C2 = I + v*C1
*
               DO J = 1, N
                  DO I = 1, M-1
                     IF( I.EQ.(M-N+J-1) ) THEN
                        C(I,J) = ONE + V(1 + (I-1)*INCV)*C(M,J)
                     ELSE
                        C(I,J) = V(1 + (I-1)*INCV)*C(M,J)
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau*v'*C2
*
*              note this is a hacky solution in place of adding
*              a {c,z}gemcv that conjugates the vector x in zgemv
*
               CALL ZGEMM('Conjugate Transpose', 'No Transpose',
     $            1, N, M-1, -TAU, V, M-1, C, LDC,
     $            ZERO, C(M,1), LDC)
*
*              C2 = C2 + v*C1
*
               CALL ZGERU(M-1, N, ONE, V, INCV, C(M,1), LDC,
     $            C, LDC)
            END IF
         END IF
      ELSE IF( RQ ) THEN
*
*        the assumed 1 in V is last in storage: ie our refector
*        is stored as [ v 1 ]
*
*        we are computing C = CH = C(I - V'TV)
*        where: V = [ v 1 ] and C = [ C2 C1 ] and T=tau is a scalar
*        with the following dimensions:
*            v\in\R^{1\times n-1}
*            C1=0\in\R^{m\times 1}
*            C2  \in\R^{m\times n-1}
*        since we are assuming that C1 is a zero column and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        this lets us simplify our above equation to get
*
*        C = CH = [ C2 0 ] - [ C2 0 ][ v' ]T[ v 1 ]
*                                    [ 1  ]
*
*        = [ C2 0 ] - C2*v'*T [ v 1 ]
*
*        = [ C2 0 ] - [ C2*v'*T*v C2*v'*T ]
*
*        = [ C2 - C2*v'*T*v, - C2*v'*T ]
*
*        so, we can order our computations as follows:
*
*        C1 = -tau*C2*v'
*        C2 = C2 + C1*v
*
*        if we also add the constraint that C2 starts as the
*        first n columns of the m-1 identity, then this simplifies into
*
*        C1 = -tau*v(n-m:n-1)'
*        C2 = I + C1*v
*
*        first check if T = tau is 0, if this is the case, all we
*        need to do is set the last row of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO I = 1, M
               C(I, N) = ZERO
            END DO
*
*           if we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF( C2I ) THEN
               CALL ZLASET('All', M, N-1, ZERO, ZERO, C, LDC)
               DO I = 1, M
                  C(I,N-M-1+I) = ONE
               END DO
            END IF
         ELSE
            IF( C2I ) THEN
*
*              C1 = -tau*v(n-m:n-1)'
*
               DO I = 1, M
                  C(I,N) = -TAU*CONJG(V(1 + (N-M-1 + I-1)*INCV))
               END DO
*
*              C2 = I + C1*v
               DO J = 1, N-1
                  DO I = 1, M
                     IF( J.EQ.(N-M + I-1) ) THEN
                        C(I,J) = ONE + C(I,N) * V(1 + (J-1)*INCV)
                     ELSE
                        C(I,J) = C(I,N) * V(1 + (J-1)*INCV)
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau*C2*v'
*
               CALL ZGEMM('No Tranpose', 'Conjugate Transpose',
     $               M, 1, N-1, -TAU, C, LDC, V, INCV,
     $               ZERO, C(1,N), LDC)
*
*              C2 = C2 + C1*v
*
               CALL ZGERU(M, N-1, ONE, C(1,N), 1, V, INCV,
     $            C, LDC)
            END IF
         END IF
      END IF
      RETURN
      END SUBROUTINE
