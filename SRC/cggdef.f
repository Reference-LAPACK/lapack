*> \brief \b CGGDEF
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> CGGDEF performs the deflation of infinite eigenvalues in the generalized
*> eigenvalue problem.
*> It applies QR and RQ factorizations to 
*> sub-blocks of A and B to properly isolate and deflate these eigenvalues.
*> \endverbatim
*
* Arguments:
* ==========
*
*> \param[in] COMPVL
*> \verbatim
*>          COMPVL is LOGICAL
*>          = .TRUE.:  Compute the left Schur vectors (VL).
*>          = .FALSE.: Do not compute the left Schur vectors.
*> \endverbatim
*>
*> \param[in] COMPVR
*> \verbatim
*>          COMPVR is LOGICAL
*>          = .TRUE.:  Compute the right Schur vectors (VR).
*>          = .FALSE.: Do not compute the right Schur vectors.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A, B, VL, and VR.
*>          N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of infinite eigenvalues to deflate (the block size).
*>          0 <= K <= N.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA, N)
*>          On entry, the matrix A to be deflated.
*>          On exit, A has been updated by the unitary transformations.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB, N)
*>          On entry, the matrix B to be deflated.
*>          On exit, B has been updated by the unitary transformations.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VL
*> \verbatim
*>          VL is COMPLEX array, dimension (LDVL, N)
*>          If COMPVL = .TRUE., the left Schur vectors are accumulated in VL.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.
*>          If COMPVL = .TRUE., LDVL >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VR
*> \verbatim
*>          VR is COMPLEX array, dimension (LDVR, N)
*>          If COMPVR = .TRUE., the right Schur vectors are accumulated in VR.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.
*>          If COMPVR = .TRUE., LDVR >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          LWORK >= max(1, 2*N).
*>          If LWORK = -1, then a workspace query is assumed;
*>          the routine
*>          only calculates the optimal size of the WORK array.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
      SUBROUTINE CGGDEF( COMPVL, COMPVR, N, K, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, WORK, LWORK, INFO )
*
* -- LAPACK-style Code --
*
* .. Scalar Arguments ..
      LOGICAL            COMPVL, COMPVR
      INTEGER            INFO, K, LDA, LDB, LDVL, LDVR, LWORK, N
* ..
* .. Array Arguments ..
      COMPLEX         A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
* ..
*
* =====================================================================
*
* .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LWKOPT
* ..
* .. External Subroutines ..
      EXTERNAL           CGEQRF, CGERQF, CUNMQR, CUNMRQ, XERBLA
* ..
* .. Intrinsic Functions ..
      INTRINSIC          REAL, CMPLX, INT, MAX
* ..
* .. Executable Statements ..
*
* Test the input parameters.
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( COMPVL .AND. LDVL.LT.MAX( 1, N ) ) THEN
         INFO = -10
      ELSE IF( COMPVR .AND. LDVR.LT.MAX( 1, N ) ) THEN
         INFO = -12
      ELSE IF( LWORK.LT.MAX( 1, 2*N ) .AND. .NOT.LQUERY ) THEN
* Minimum workspace is N elements for TAU arrays + N elements for routines
         INFO = -14
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
* Compute optimal workspace
*
         LWKOPT = MAX( 1, 2*N )
         IF( K.GT.0 ) THEN
            CALL CGEQRF( N, K, A, LDA, WORK, WORK, -1, IINFO )
            LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
*
            IF( N.GT.K ) THEN
               CALL CUNMQR( 'L', 'C', N, N - K, K, A, LDA, WORK,
     $                      B( 1, K + 1 ), LDB, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
            END IF
*
            IF( COMPVL ) THEN
               CALL CUNMQR( 'R', 'N', N, N, K, A, LDA, WORK, VL,
     $                      LDVL, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
            END IF
         END IF
*
         IF( N.GT.K ) THEN
            CALL CGERQF( N - K, N - K, B( K + 1, K + 1 ), LDB, WORK,
     $                   WORK, -1, IINFO )
            LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
*
            IINFO = -1
            CALL CUNMRQ( 'R', 'C', N, N - K, N - K, 
     $                   B( K + 1, K + 1 ),
     $                   LDB, WORK, A( 1, K + 1 ), LDA, WORK, -1,
     $                   IINFO )
            LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
*
            IF( K.GT.0 ) THEN
               CALL CUNMRQ( 'R', 'C', K, N - K, N - K,
     $                      B( K + 1, K + 1 ), LDB, WORK,
     $                      B( 1, K + 1 ), LDB, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
            END IF
*
            IF( COMPVR ) THEN
               CALL CUNMRQ( 'R', 'C', N, N - K, N - K,
     $                      B( K + 1, K + 1 ), LDB, WORK,
     $                      VR( 1, K + 1 ), LDVR, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( REAL( WORK( 1 ) ) ) )
            END IF
         END IF
*
         WORK( 1 ) = CMPLX( REAL( LWKOPT ), 0.0E0 )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CGGDEF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
* Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
* QR portion
* Factorize A(1:N, 1:K)
* WORK(1:K) holds the TAU scalar factors for A
*
      IF( K.GT.0 ) THEN
         CALL CGEQRF( N, K, A, LDA, WORK( 1 ), WORK( N + 1 ),
     $                LWORK - N, IINFO )
*
* Apply Q^H to A(1:N, K+1:N)
*
         IF( N.GT.K ) THEN
            CALL CUNMQR( 'L', 'C', N, N - K, K, A, LDA, WORK( 1 ),
     $                   A( 1, K + 1 ), LDA, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
*
* Apply Q^H to B(1:N, K+1:N)
*
            CALL CUNMQR( 'L', 'C', N, N - K, K, A, LDA, WORK( 1 ),
     $                   B( 1, K + 1 ), LDB, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
         END IF
*
         IF( COMPVL ) THEN
* Accumulate Q in VL: VL = VL * Q
            CALL CUNMQR( 'R', 'N', N, N, K, A, LDA, WORK( 1 ), VL,
     $                   LDVL, WORK( N + 1 ), LWORK - N, IINFO )
         END IF
      END IF
*
* RQ portion
* Factorize B(K+1:N, K+1:N)
* WORK(K+1:N) holds the TAU scalar factors for B
*
      IF( N.GT.K ) THEN
         CALL CGERQF( N - K, N - K, B( K + 1, K + 1 ), LDB,
     $                WORK( K + 1 ), WORK( N + 1 ), LWORK - N, IINFO )
*
* Apply Q^H to A(1:N, K+1:N) from the right: A = A * Q^H
*
         CALL CUNMRQ( 'R', 'C', N, N - K, N - K, B( K + 1, K + 1 ),
     $                LDB, WORK( K + 1 ), A( 1, K + 1 ), LDA,
     $                WORK( N + 1 ), LWORK - N, IINFO )
*
* Apply Q^H to B(1:K, K+1:N) from the right: B = B * Q^H
*
         IF( K.GT.0 ) THEN
            CALL CUNMRQ( 'R', 'C', K, N - K, N - K,
     $                   B( K + 1, K + 1 ), LDB, WORK( K + 1 ),
     $                   B( 1, K + 1 ), LDB, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
         END IF
*
         IF( COMPVR ) THEN
* Accumulate Q^H in VR: VR(1:N, K+1:N) = VR(1:N, K+1:N) * Q^H
            CALL CUNMRQ( 'R', 'C', N, N - K, N - K,
     $                   B( K + 1, K + 1 ), LDB, WORK( K + 1 ),
     $                   VR( 1, K + 1 ), LDVR, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
         END IF
      END IF
*
* Zero out strictly lower triangular parts
*
      DO J = 1, K
         DO I = J + 1, N
            A( I, J ) = ( 0.0E0, 0.0E0 )
         END DO
      END DO
*
      DO J = K + 1, N
         DO I = J + 1, N
            B( I, J ) = ( 0.0E0, 0.0E0 )
         END DO
      END DO
*
      RETURN
*
* End of CGGDEF
*
      END