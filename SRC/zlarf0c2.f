*> \brief \b ZLARF0C2 applies an elementary reflector to a rectangular matrix
*> with a 0 row/column while constructing the explicit Q factor.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*     SUBROUTINE ZLARF0C2( C2JOB, SIDE, DIRECT, STOREV, M, N, TAU,
*                          V, C, LDC )
*
*     .. Scalar Arguments ..
*     INTEGER           M, N, LDC
*     CHARACTER         C2JOB, SIDE, DIRECT, STOREV
*     COMPLEX*16        TAU
*     ..
*     .. Array Arguments ..
*     COMPLEX*16        C(LDC,*), V( * )
*     ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORG1R generates the matrix associated with a householder reflector,
*> which is defined as the first n columns of a single elementary
*> reflector of order m
*>
*>       H = I - \tau*v*v**H
*>
*> as returned by ZLARFG.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] C2JOB
*> \verbatim
*>          C2JOB is CHARACTER*1
*>          = 'I': Assume C2 is the identity matrix
*>          Otherwise: Treat C2 as a general matrix that we reference
*> \endverbatim
*>
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'L': apply H or H**T from the Left
*>          = 'R': apply H or H**T from the Right
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
*>          The number of rows of the matrix H. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix H. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          TAU must contain the value TAU returned from ZLARFG
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is COMPLEX*16 array, dimension
*>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*>          The vector v in the representation of H. V is not used if
*>          TAU = 0. See Further Details.
*> \endverbatim
*>
*> \param[in] INCV
*> \verbatim
*>          INCV is INTEGER
*>          The increment between elements of v. INCV <> 0.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is COMPLEX*16 array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by H*C or C*H.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument has an illegal value
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
*>  The storage of the vector which define H is best illustrated with
*>  the following example. We use a v of length 5, (this is M if
*>  STOREV = 'C' and N if STOREV = 'R')
*>
*>  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*>
*>               V = (  1 )                 V = (  1 v2 v3 v4 v5 )
*>                   ( v2 )
*>                   ( v3 )
*>                   ( v4 )
*>                   ( v5 )
*>
*>  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*>
*>               V = ( v1 )                 V = ( v1 v2 v3 v4  1 )
*>                   ( v2 )
*>                   ( v3 )
*>                   ( v4 )
*>                   (  1 )
*>
*>  Also see zlarf1f and zlarf1l
*>  for a similar routine but without the assumed 0 block in C
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE ZLARF0C2( C2JOB, SIDE, DIRECT, STOREV, M, N, TAU,
     $                     V, INCV, C, LDC )
*
*     .. Scalar Arguments ..
      INTEGER           M, N, LDC, INCV
      CHARACTER         C2JOB, SIDE, DIRECT, STOREV
      COMPLEX*16        TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16        C(LDC,*), V( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16        ONE, ZERO
      PARAMETER         ( ONE = (1.0D+0, 0.0D+0),
     $                    ZERO = (0.0D+0, 0.0D+0) )
*     ..
*     .. Local Scalars ..
      LOGICAL           QR, LQ, QL, RQ, SIDEL, SIDER, COLV, C2I, DIRF
      INTEGER           I, J
*     ..
*     .. External Functions ..
      LOGICAL           LSAME
      EXTERNAL          LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL          ZGEMV, ZGERU, ZLASET, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         CONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible. When \tau is 0, H is defined to be I,
*     so the application doesn't change C at all. However we still flow
*     through as we will need to set the appropriate row of C to the 0 row
*     In addition, our blas kernels should do nothing if the input scalar
*     is 0.
*
      IF( M.LE.0.OR.N.LE.0 ) THEN
         RETURN
      END IF
      DIRF  = LSAME(DIRECT, 'F')
      COLV  = LSAME(STOREV, 'C')
      C2I   = LSAME(C2JOB,  'I')
      SIDEL = LSAME(SIDE,   'L')
      SIDER = LSAME(SIDE,   'R')
      IF( .NOT.(SIDEL.OR.SIDER) ) THEN
         CALL XERBLA('ZLARF0C2', 2) ! For consistency with dlarfb0c2
      END IF
*
*     Determine which of the 4 modes we are doing.
*     QR is when we have the first element of v being 1, and is stored
*     as a column vector
*
      QR = DIRF.AND.COLV
*
*     LQ is when we have the first element of v being 1, and is stored
*     as a row vector
*
      LQ = DIRF.AND.(.NOT.COLV)
*
*     QL is when we have the last element of v being 1, and is stored
*     as a column vector
*
      QL = (.NOT.DIRF).AND.COLV
*
*     RQ is when we have the last element of v being 1, and is stored
*     as a row vector
*
      RQ = (.NOT.DIRF).AND.(.NOT.COLV)
      IF( QR ) THEN
*        We are computing C = HC = (I - VTV')C
*        Where: V = [ 1 ], C = [ C1 ], and T=tau is a scalar
*                   [ v ]      [ C2 ]
*        with the following dimensions:
*            v\in\R^{m-1\times 1}
*            C1=0\in\R^{1\times n}
*            C2  \in\R^{m-1\times n}
*        Since we are assuming that C1 is a zero row and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        This lets us simplify our above equation to get
*
*        C = HC = (I - [  1 ]T [ 1,  v'])[ 0  ]
*                      [  v ]            [ C2 ]
*          = [ 0  ] - [  1 ]T* v'*C2
*            [ C2 ]   [  v ]
*
*          = [ 0  ] - [ T* v'*C2 ]
*            [ C2 ]   [  v*T* v'*C2 ]
*
*          = [        -T* v'*C2 ]
*            [ C2 -  v*T* v'*C2 ]
*
*        So, we can order our computations as follows:
*
*        C1 = -tau * v'*C2
*        C2 = C2 + v*C1
*
*        If we also add the constraint that C2 starts as the
*        first n columns of the m-1 identity, then this simplifies into
*
*        C1 = -tau*v(1:n)'
*        C2 = I + v*C1
*
*        First check if T = tau is 0, if this is the case, all we
*        need to do is set the first row of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO J = 1, N
               C(1,J) = ZERO
            END DO
*
*           If we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF ( C2I ) THEN
               CALL ZLASET('All', M-1, N, ZERO, ONE, C(2,1), LDC)
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
*              C1 = v'*C2 = -tau*(C2**T*conj(v))**T
*
               CALL ZGEMCV('Transpose', M-1, N, -TAU,
     $              C(2,1), LDC, V, INCV, ZERO, C, LDC)
*
*              C2 = C2 + v*C1
*
               CALL ZGERU(M-1, N, ONE, V, INCV, C, LDC, C(2,1), LDC)
            END IF
         END IF
      ELSE IF( LQ ) THEN
*
*        We are computing C = CH = C(I-V'T'V)
*        Where: V = [  1 v ], C = [ 0 C2 ], and T=tau is a scalar
*        with the following dimensions:
*            v\in\R^{1\times n-1}
*            C1=0\in\R^{m\times 1}
*            C2  \in\R^{m\times n-1}
*        Since we are assuming that C1 is a zero row and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        This lets us simplify our above equation to get
*
*        C = CH = [ 0 C2 ] - [ 0 C2 ] [ 1  ] T [ 1 v ]
*                                     [ v' ]
*
*        = [ 0 C2 ] - C2*v'*T*[ 1 v ]
*
*        = [ 0 C2 ] - [ C2*v'*T C2*v'*T*v ]
*
*        = [ -C2*v'*T   C2 - C2*v'*T*v ]
*
*        This means we can order our computation as follows:
*
*        C1 = -tau*C2*v'
*        C2 = C2 + C1*v
*
*        If we also add the constraint that C2 starts as the
*        first m columns of the n-1 identity, then this simplifies into
*
*        C1 = -tau*v(1:m)'
*        C2 = I + C1*v
*
*        First check if T = tau is 0, if this is the case, all we
*        need to do is set the first column of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO I = 1, M
               C(I,1) = ZERO
            END DO
*
*           If we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF ( C2I ) THEN
               CALL ZLASET('All', M, N-1, ZERO, ONE, C(1,2), LDC)
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
*              no routines that perform this operation exist, so
*              we compute this columnwise
*
               DO J = 1, N-1
                  DO I = 1, M
                     IF( I.EQ.J ) THEN
                        C(I,1+J) = ONE + C(I,1)*V(1 + (J-1)*INCV)
                     ELSE
                        C(I,1+J) = C(I,1)*V(1 + (J-1)*INCV)
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau*C2*v'
*
               CALL ZGEMCV('No Transpose', M, N-1, -TAU,
     $            C(1,2), LDC, V, INCV, ZERO, C, 1)
*
*              C2 = C2 + C1*v
*
               CALL ZGERU(M, N-1, ONE, C, 1, V, INCV, C(1,2), LDC)
            END IF
         END IF
      ELSE IF( QL ) THEN
*
*        We are computing C = HC = (I - VTV')C
*        Where: V = [ v ] and C = [ C2 ] and T=tau is a scalar
*                   [ 1 ]         [ C1 ]
*        with the following dimensions:
*            v\in\R^{m-1\times 1}
*            C1=0\in\R^{1\times n}
*            C2  \in\R^{m-1\times n}
*        Since we are assuming that C1 is a zero row and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        This lets us simplify our above equation to get
*
*        C = HC = [ C2 ] - [ v ]T[ v' 1 ][ C2 ]
*                 [ 0  ]   [ 1 ]         [ 0  ]
*
*        = [ C2 ] - [ v ]T*v'C2
*          [ 0  ]   [ 1 ]
*
*        = [ C2 - v*T*v'*C2 ]
*          [    - T*v'*C2 ]
*
*        So, we can order our computations as follows:
*
*        C1 = -tau*C2**T*conjg(v)
*        C2 = C2 + v*C1
*
*        If we also add the constraint that C2 starts as the
*        first n columns of the m-1 identity, then this simplifies into
*
*        C1 = -tau*conjg(v(m-n:m-1))
*        C2 = I + v*C1
*
*        First check if T = tau is 0, if this is the case, all we
*        need to do is set the last row of C to be 0 then exit
*
         IF( TAU.EQ.ZERO ) THEN
            DO J = 1, N
               C(M, J) = ZERO
            END DO
*
*           If we assumed C2 was the identity matrix, we must explicitly
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
*              C1 = -tau*conjg(v(m-n:m-1))
*
               DO J = 1, N
                  C(M,J) = -TAU*CONJG(V(1 + (M-1-N+J-1)*INCV))
               END DO
*
*              C2 = I + v*C1
*
               DO J = 1, N
                  DO I = 1, M-1
                     IF( I.EQ.(M-N-1+J) ) THEN
                        C(I,J) = ONE + V(1 + (I-1)*INCV)*C(M,J)
                     ELSE
                        C(I,J) = V(1 + (I-1)*INCV)*C(M,J)
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau*C2**T*conjg(v)
*
               CALL ZGEMCV('Transpose', M-1, N, -TAU, C, LDC,
     $               V, INCV, ZERO, C(M,1), LDC)
*
*              C2 = C2 + v*C1
*
               CALL ZGERU(M-1, N, ONE, V, INCV, C(M,1), LDC, C, LDC)
            END IF
         END IF
      ELSE IF( RQ ) THEN
*
*        We are computing C = CH = (I - V'TV)C
*        Where: V = [ v 1 ] and C = [ C2 C1 ] and T=tau is a scalar
*        with the following dimensions:
*            v\in\R^{1\times n-1}
*            C1=0\in\R^{m\times 1}
*            C2  \in\R^{m\times n-1}
*        Since we are assuming that C1 is a zero column and it will be
*        overwritten on exit, we can use this spot as a temporary workspace
*        without having to allocate anything extra.
*        This lets us simplify our above equation to get
*
*        C = CH = [ C2, 0 ] - [ C2, 0 ][ v' ]T[ v 1 ]
*                                      [ 1  ]
*
*        = [ C2, 0 ] - C2v'T[ v 1 ]
*
*        = [ C2 - C2v'Tv, -C2v'T ]
*
*        So, we can order our computations as follows:
*
*        C1 = -tau*C2v'
*        C2 = C2 + C1*v
*
*        If we also add the constraint that C2 starts as the first
*        m rows of the n-1 identity, then this simplifies into
*
*        C1 = -tau*conjg(v(n-m:n-1))
*        C2 = C2 + C1*v
*
         IF( TAU.EQ.ZERO ) THEN
            DO I = 1, M
               C(I, N) = ZERO
            END DO
*
*           If we assumed C2 was the identity matrix, we must explicitly
*           set this before we exit
*
            IF( C2I ) THEN
               CALL ZLASET('All', M, N-1, ZERO, ZERO, C, LDC)
               DO I = 1, M
                  C(I, N-M-1+I) = ONE
               END DO
            END IF
         ELSE
            IF( C2I ) THEN
*
*              C1 = -tau*v(n-m:n-1)
*
               DO I = 1, M
                  C(I,N) = -TAU*CONJG(V(1 + (N-1-M+I-1)*INCV))
               END DO
*
*              C2 = I + C1*V
*
               DO J = 1, N-1
                  DO I = 1, M
                     IF( J.EQ.(N-M-1+I) ) THEN
                        C(I,J) = ONE + C(I, N)*V(1 + (J-1)*INCV)
                     ELSE
                        C(I,J) = C(I, N)*V(1 + (J-1)*INCV)
                     END IF
                  END DO
               END DO
            ELSE
*
*              C1 = -tau*C2v'
*
               CALL ZGEMCV('No Transpose', M, N-1, -TAU, C, LDC,
     $            V, INCV, ZERO, C(1,N), 1)
*
*              C2 = C2 + C1*v
*
               CALL ZGERU(M, N-1, ONE, C(1,N), 1, V, INCV, C, LDC)
            END IF
         END IF
      END IF
      RETURN
      END SUBROUTINE
