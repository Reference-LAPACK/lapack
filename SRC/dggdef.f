*> \brief \b DGGDEF
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> DGGDEF performs the deflation of infinite eigenvalues in the generalized
*> eigenvalue problem. It applies QR and RQ factorizations to 
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
*>          The order of the matrices A, B, VL, and VR.  N >= 0.
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
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the matrix A to be deflated.
*>          On exit, A has been updated by the orthogonal transformations.
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
*>          B is DOUBLE PRECISION array, dimension (LDB, N)
*>          On entry, the matrix B to be deflated.
*>          On exit, B has been updated by the orthogonal transformations.
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
*>          VL is DOUBLE PRECISION array, dimension (LDVL, N)
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
*>          VR is DOUBLE PRECISION array, dimension (LDVR, N)
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
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1, 2*N).
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
      SUBROUTINE DGGDEF( COMPVL, COMPVR, N, K, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, WORK, LWORK, INFO )
*
* -- LAPACK-style Code --
*
* .. Scalar Arguments ..
      LOGICAL            COMPVL, COMPVR
      INTEGER            INFO, K, LDA, LDB, LDVL, LDVR, LWORK, N
* ..
* .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
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
      EXTERNAL           DGEQRF, DGERQF, DORMQR, DORMRQ, XERBLA
* ..
* .. Intrinsic Functions ..
      INTRINSIC          DBLE, INT, MAX
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
            CALL DGEQRF( N, K, A, LDA, WORK, WORK, -1, IINFO )
            LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
*
            IF( N.GT.K ) THEN
               CALL DORMQR( 'L', 'T', N, N - K, K, A, LDA, WORK,
     $                      B( 1, K + 1 ), LDB, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
            END IF
*
            IF( COMPVL ) THEN
               CALL DORMQR( 'R', 'N', N, N, K, A, LDA, WORK, VL,
     $                      LDVL, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
            END IF
         END IF
*
         IF( N.GT.K ) THEN
            CALL DGERQF( N - K, N - K, B( K + 1, K + 1 ), LDB, WORK,
     $                   WORK, -1, IINFO )
            LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
*
            IINFO = -1
            CALL DORMRQ( 'R', 'T', N, N - K, N - K, 
     $                   B( K + 1, K + 1 ),
     $                   LDB, WORK, A( 1, K + 1 ), LDA, WORK, -1,
     $                   IINFO )
            LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
*
            IF( K.GT.0 ) THEN
               CALL DORMRQ( 'R', 'T', K, N - K, N - K,
     $                      B( K + 1, K + 1 ), LDB, WORK,
     $                      B( 1, K + 1 ), LDB, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
            END IF
*
            IF( COMPVR ) THEN
               CALL DORMRQ( 'R', 'T', N, N - K, N - K,
     $                      B( K + 1, K + 1 ), LDB, WORK,
     $                      VR( 1, K + 1 ), LDVR, WORK, -1, IINFO )
               LWKOPT = MAX( LWKOPT, N + INT( WORK( 1 ) ) )
            END IF
         END IF
*
         WORK( 1 ) = DBLE( LWKOPT )
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGDEF', -INFO )
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
         CALL DGEQRF( N, K, A, LDA, WORK( 1 ), WORK( N + 1 ),
     $                LWORK - N, IINFO )
*
* Apply Q^T to A(1:N, K+1:N)
*
         IF( N.GT.K ) THEN
            CALL DORMQR( 'L', 'T', N, N - K, K, A, LDA, WORK( 1 ),
     $                   A( 1, K + 1 ), LDA, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
*
* Apply Q^T to B(1:N, K+1:N)
*
            CALL DORMQR( 'L', 'T', N, N - K, K, A, LDA, WORK( 1 ),
     $                   B( 1, K + 1 ), LDB, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
         END IF
*
         IF( COMPVL ) THEN
* Accumulate Q in VL: VL = VL * Q
            CALL DORMQR( 'R', 'N', N, N, K, A, LDA, WORK( 1 ), VL,
     $                   LDVL, WORK( N + 1 ), LWORK - N, IINFO )
         END IF
      END IF
*
* RQ portion
* Factorize B(K+1:N, K+1:N)
* WORK(K+1:N) holds the TAU scalar factors for B
*
      IF( N.GT.K ) THEN
         CALL DGERQF( N - K, N - K, B( K + 1, K + 1 ), LDB,
     $                WORK( K + 1 ), WORK( N + 1 ), LWORK - N, IINFO )
*
* Apply Q^T to A(1:N, K+1:N) from the right: A = A * Q^T
*
         CALL DORMRQ( 'R', 'T', N, N - K, N - K, B( K + 1, K + 1 ),
     $                LDB, WORK( K + 1 ), A( 1, K + 1 ), LDA,
     $                WORK( N + 1 ), LWORK - N, IINFO )
*
* Apply Q^T to B(1:K, K+1:N) from the right: B = B * Q^T
*
         IF( K.GT.0 ) THEN
            CALL DORMRQ( 'R', 'T', K, N - K, N - K,
     $                   B( K + 1, K + 1 ), LDB, WORK( K + 1 ),
     $                   B( 1, K + 1 ), LDB, WORK( N + 1 ),
     $                   LWORK - N, IINFO )
         END IF
*
         IF( COMPVR ) THEN
* Accumulate Q^T in VR: VR(1:N, K+1:N) = VR(1:N, K+1:N) * Q^T
            CALL DORMRQ( 'R', 'T', N, N - K, N - K,
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
            A( I, J ) = 0.0D0
         END DO
      END DO
*
      DO J = K + 1, N
         DO I = J + 1, N
            B( I, J ) = 0.0D0
         END DO
      END DO
*
      RETURN
*
* End of DGGDEF
*
      END