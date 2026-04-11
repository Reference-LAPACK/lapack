*> \brief \b DGGHD4
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> DGGHD4 reduces a pair of real matrices (A,B) to generalized upper
*> Hessenberg form using orthogonal transformations, where A is a
*> Hessenberg matrix and B is upper triangular. This routine
*> employs a blocked algorithm (Steel et al) to improve efficiency
*> and scalability.
*> \endverbatim
*
* Arguments:
* ==========
*
*> \param[in] COMPQ
*> \verbatim
*>          COMPQ is CHARACTER*1
*>          = 'N': do not compute Q;
*>          = 'I': Q is initialized to the unit matrix, and the
*>                 orthogonal matrix Q is returned;
*>          = 'V': Q must contain an orthogonal matrix Q1 on entry,
*>                 and the product Q1*Q is returned.
*> \endverbatim
*>
*> \param[in] COMPZ
*> \verbatim
*>          COMPZ is CHARACTER*1
*>          = 'N': do not compute Z;
*>          = 'I': Z is initialized to the unit matrix, and the
*>                 orthogonal matrix Z is returned;
*>          = 'V': Z must contain an orthogonal matrix Z1 on entry,
*>                 and the product Z1*Z is returned.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in] ILO
*> \param[in] IHI
*> \verbatim
*>          ILO and IHI are INTEGER
*>          It is assumed that A is already upper triangular in rows
*>          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*>          set by a previous call to DGGBAL.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA, N)
*>          On entry, the N-by-N general matrix to be reduced.
*>          On exit, the upper triangle and the first subdiagonal of A
*>          are overwritten with the upper Hessenberg matrix H.
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
*>          On entry, the N-by-N upper triangular matrix B.
*>          On exit, the upper triangular matrix T = Q**T B Z.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] Q
*> \verbatim
*>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
*>          If COMPQ='N', then Q is not referenced.
*>          If COMPQ='V', on entry Q must contain an orthogonal matrix.
*>          If COMPQ='I', on exit Q contains the orthogonal matrix.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.
*>          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*> \endverbatim
*>
*> \param[in,out] Z
*> \verbatim
*>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
*>          If COMPZ='N', then Z is not referenced.
*>          If COMPZ='V', on entry Z must contain an orthogonal matrix.
*>          If COMPZ='I', on exit Z contains the orthogonal matrix.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.
*>          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
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
*>          The dimension of the array WORK.
*>          If LWORK = -1, then a workspace query is assumed.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim
      SUBROUTINE DGGHD4( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB,
     $                   Q, LDQ, Z, LDZ, WORK, LWORK, INFO )
*
* -- LAPACK-like routine --
*
* .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
* ..
* .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   WORK( * ), Z( LDZ, * )
* ..
*
* =====================================================================
*
* .. Local Scalars ..
      LOGICAL            ILQ, ILZ, LQUERY
      INTEGER            I, J, IERR, IWORK, K, LWKOPT, NCOLS,
     $                   NITER, COL_LEN, LMIN
      DOUBLE PRECISION   EPS, NORM, NORMA, TOL, AMAX
* ..
* .. Local Arrays ..
      DOUBLE PRECISION   DUMMY( 1 )
* ..
* .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANGE
      EXTERNAL           LSAME, DLAMCH, DLANGE
* ..
* .. External Subroutines ..
      EXTERNAL           DGEHRD, DGERQF, DLACPY, DLASET, DORMHR,
     $                   DORMRQ, DTRSM, XERBLA
* ..
* .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN

      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
* ..
* .. Executable Statements ..
*
* Decode and test the input parameters.
*
      INFO = 0
      ILQ = LSAME( COMPQ, 'V' ) .OR. LSAME( COMPQ, 'I' )
      ILZ = LSAME( COMPZ, 'V' ) .OR. LSAME( COMPZ, 'I' )
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.LSAME( COMPQ, 'N' ) .AND. .NOT.ILQ ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( COMPZ, 'N' ) .AND. .NOT.ILZ ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 .OR. ILO.GT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( IHI.LT.MIN( ILO, N ) .OR. IHI.GT.N ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ILQ .AND. LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -11
      ELSE IF( ILZ .AND. LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -13
      END IF
*
* Compute optimal workspace.
*
      IF( INFO.EQ.0 ) THEN
         LWKOPT = 2*( N + 1 )
         CALL DGEHRD( N, 1, N, A, LDA, DUMMY, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         CALL DORMHR( 'L', 'T', N, N, 1, N, A, LDA, DUMMY, A, LDA,
     $                WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         CALL DGERQF( N, N, A, LDA, DUMMY, WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         CALL DORMRQ( 'R', 'T', N, N, N, A, LDA, DUMMY, A, LDA,
     $                WORK, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( WORK( 1 ) ) )
         LWKOPT = LWKOPT + N*( N + 1 )
         LMIN = LWKOPT
         WORK( 1 ) = DBLE( LWKOPT )
         IF( LWORK.LT.LMIN .AND. .NOT.LQUERY ) THEN
            INFO = -15
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGGHD4', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
* Initialize Q and Z if requested.
*
      IF( LSAME( COMPQ, 'I' ) ) THEN
         CALL DLASET( 'F', N, N, 0.0D0, 1.0D0, Q, LDQ )
      END IF
      IF( LSAME( COMPZ, 'I' ) ) THEN
         CALL DLASET( 'F', N, N, 0.0D0, 1.0D0, Z, LDZ )
      END IF
*
* Quick return if possible
*
      IF( N.LE.1 ) RETURN
*
      NORMA = DLANGE( 'F', N, N, A, LDA, WORK )
      EPS = DLAMCH( 'E' )
      TOL = EPS * NORMA
      NITER = 0
      K = ILO
      IWORK = N*N + 1

      AMAX = ZERO
      DO J = 1, N
         DO I = 1, N
            AMAX = MAX( AMAX, ABS(A(I, J)) )
         END DO
      END DO
      IF (AMAX.EQ.ZERO) THEN
         AMAX = ONE
      END IF
      TOL = EPS * AMAX

*
* Outer reduction loop: iterate until K reaches IHI or NITER hits 50.
*
      DO WHILE( K.LT.IHI .AND. NITER.LT.50 )
*
* Inner deflation loop: advance K while subdiagonal entries are small.
*
         DO WHILE( K.LT.IHI )
            COL_LEN = IHI - K - 1
            
            IF( COL_LEN.LE.0 ) THEN
               K = K + 1
               CYCLE
            END IF
            
            NORM = 0.0D0
            DO I = 1, COL_LEN
               NORM = MAX( NORM, ABS( A( K + 1 + I, K ) ) )
            END DO
            
            IF( NORM.GT.TOL ) THEN
               EXIT
            END IF

*
* Deflate strictly bounded values
*
            DO I = 1, COL_LEN
               A( K + 1 + I, K ) = 0.0D0
            END DO
            K = K + 1
         END DO

         IF (K.GE.IHI) THEN
            EXIT
         END IF

         NITER = NITER + 1
         NCOLS = IHI - K + 1
*
* Copy A(K:IHI, K:IHI) into X (stored at WORK(1))
*
         CALL DLACPY( 'A', NCOLS, NCOLS, A( K, K ), LDA, WORK, NCOLS )
*
* DTRSM to solve X * B = A
*
         CALL DTRSM( 'R', 'U', 'N', 'N', NCOLS, NCOLS, ONE / AMAX,
     $               B( K, K ), LDB, WORK, NCOLS )
*
* Hessenberg reduction of X
*
         CALL DGEHRD( NCOLS, 1, NCOLS, WORK, NCOLS, WORK( IWORK ),
     $                WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1, IERR )
*
* Apply Qc**T to A and B from the left
*
         CALL DORMHR( 'L', 'T', NCOLS, NCOLS, 1, NCOLS, WORK, NCOLS,
     $                WORK( IWORK ), A( K, K ), LDA,
     $                WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1, IERR )
         CALL DORMHR( 'L', 'T', NCOLS, NCOLS, 1, NCOLS, WORK, NCOLS,
     $                WORK( IWORK ), B( K, K ), LDB,
     $                WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1, IERR )
         IF( ILQ ) THEN
            CALL DORMHR( 'R', 'N', N, NCOLS, 1, NCOLS, WORK, NCOLS,
     $                   WORK( IWORK ), Q( 1, K ), LDQ,
     $                   WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1,
     $                   IERR )
         END IF

*
* RQ factorization of B(K:IHI, K:IHI)
*
         CALL DGERQF( NCOLS, NCOLS, B( K, K ), LDB, WORK( IWORK ),
     $                WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1, IERR )
*
* Apply Zc**T to A from the right (Updates rows 1:IHI)
*
         CALL DORMRQ( 'R', 'T', IHI, NCOLS, NCOLS, B( K, K ), LDB,
     $                WORK( IWORK ), A( 1, K ), LDA,
     $                WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1, IERR )
*
* Apply Zc**T to B from the right (CRITICAL FIX: Updates rows 1:K-1)
*
         CALL DORMRQ( 'R', 'T', K-1, NCOLS, NCOLS, B( K, K ), LDB,
     $                WORK( IWORK ), B( 1, K ), LDB,
     $                WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1, IERR )
         IF( ILZ ) THEN
            CALL DORMRQ( 'R', 'T', N, NCOLS, NCOLS, B( K, K ), LDB,
     $                   WORK( IWORK ), Z( 1, K ), LDZ,
     $                   WORK( IWORK+NCOLS ), LWORK-IWORK-NCOLS+1,
     $                   IERR )
         END IF
*
* Zero lower triangle of B globally
*
         CALL DLASET( 'L', N-1, N-1, 0.0D0, 0.0D0, B( 2, 1 ), LDB )

*
      END DO
*
* End of DGGHD4
*
      END SUBROUTINE DGGHD4