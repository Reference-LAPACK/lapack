*> \brief \b SGGPRP
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> SGGPRP is a preprocessing routine for the generalized 
*> eigenvalue problem. It prepares the matrices A and B by applying 
*> an initial Rank-Revealing QR (RRQR) factorization on B, transposing
*> and reversing arrays, and initializing orthogonal transformation 
*> matrices Q and Z if requested.
*> \endverbatim
*
* Arguments:
* ==========
*
*> \param[in] WANTQ
*> \verbatim
*>          WANTQ is LOGICAL
*>          = .TRUE.:  Initialize and compute the orthogonal matrix Q.
*>          = .FALSE.: Do not compute Q.
*> \endverbatim
*>
*> \param[in] WANTZ
*> \verbatim
*>          WANTZ is LOGICAL
*>          = .TRUE.:  Initialize and compute the orthogonal matrix Z.
*>          = .FALSE.: Do not compute Z.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the matrix A.
*>          On exit, A is overwritten by the preprocessed matrix.
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
*>          B is REAL array, dimension (LDB, N)
*>          On entry, the matrix B.
*>          On exit, B is overwritten by the preprocessed upper 
*>          triangular matrix from the RRQR factorization.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] Q
*> \verbatim
*>          Q is REAL array, dimension (LDQ, N)
*>          If WANTQ = .TRUE., contains the initialized orthogonal matrix.
*> \endverbatim
*>
*> \param[in] LDQ
*> \verbatim
*>          LDQ is INTEGER
*>          The leading dimension of the array Q.
*>          If WANTQ = .TRUE., LDQ >= max(1,N).
*> \endverbatim
*>
*> \param[out] Z
*> \verbatim
*>          Z is REAL array, dimension (LDZ, N)
*>          If WANTZ = .TRUE., contains the initialized orthogonal matrix.
*> \endverbatim
*>
*> \param[in] LDZ
*> \verbatim
*>          LDZ is INTEGER
*>          The leading dimension of the array Z.
*>          If WANTZ = .TRUE., LDZ >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the norm of B, and 
*>          WORK(2) returns the tolerance computed.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK.
*>          If LWORK = -1, then a workspace query is assumed.
*> \endverbatim
*>
*> \param[out] JPVT
*> \verbatim
*>          JPVT is INTEGER array, dimension (N)
*>          On exit, contains the pivot indices from SGEQP3.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit.
*>          < 0:  if INFO = -i, the i-th argument had an illegal value.
*> \endverbatim      
      SUBROUTINE SGGPRP( WANTQ, WANTZ, N, A, LDA, B, LDB, Q, LDQ,
     $                   Z, LDZ, WORK, LWORK, JPVT, INFO )
*
* -- LAPACK-style preprocessing routine --
*
* .. Scalar Arguments ..
      LOGICAL            WANTQ, WANTZ
      INTEGER            INFO, LDA, LDB, LDQ, LDZ, LWORK, N
* ..
* .. Array Arguments ..
      INTEGER            JPVT( * )
      REAL   A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   WORK( * ), Z( LDZ, * )
* ..
*
* =====================================================================
*
* .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IERR, J, LWKOPT, WORKNEEDED
      REAL   NORMB, TOL
* ..
* .. Local Arrays ..
      REAL   DUMMY( 1 )
* ..
* .. External Subroutines ..
      EXTERNAL           SGEQP3, SLAREV, SLAPMR, SLATRN, SORGQR, SORMQR,
     $                   XERBLA
* ..
* .. External Functions ..
      REAL   SLAMCH, SLANGE
      EXTERNAL           SLAMCH, SLANGE
* ..
* .. Intrinsic Functions ..
      INTRINSIC          REAL, INT, MAX
* ..
* .. Executable Statements ..
*
* Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( WANTQ .AND. LDQ.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( WANTZ .AND. LDZ.LT.MAX( 1, N ) ) THEN
         INFO = -11
      END IF
*
* Compute workspace if no argument errors occurred
*
      IF( INFO.EQ.0 ) THEN
         LWKOPT = 0
         CALL SGEQP3( N, N, B, LDB, JPVT, DUMMY, DUMMY, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( DUMMY( 1 ) ) )
*
         CALL SORMQR( 'R', 'N', N, N, N, B, LDB, DUMMY, A, LDA,
     $                DUMMY, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( DUMMY( 1 ) ) )
*
         CALL SORGQR( N, N, N, Z, LDZ, DUMMY, DUMMY, -1, IERR )
         LWKOPT = MAX( LWKOPT, INT( DUMMY( 1 ) ) )
*
* Total workspace must accommodate TAU (size N) plus the optimal
* workspace sizes requested by the underlying routines.
*
         WORKNEEDED = N + LWKOPT
*
         IF( LWORK.LT.WORKNEEDED .AND. .NOT.LQUERY ) THEN
            INFO = -13
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SGGPRP', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         WORK( 1 ) = REAL( WORKNEEDED )
         RETURN
      END IF
*
* Compute norm of B and tolerance
*
      NORMB = SLANGE( 'F', N, N, B, LDB, WORK )
      TOL = NORMB * SLAMCH( 'E' )
*
* B = B^T
*
      CALL SLATRN( N, B, LDB )
*
* B = B[:, ::-1]
*
      CALL SLAREV( 'C', N, N, B, LDB )
*
* RRQR(B)
*
      DO I = 1, N
         JPVT( I ) = 0
      END DO
      CALL SGEQP3( N, N, B, LDB, JPVT, WORK, WORK( N + 1 ),
     $             LWORK - N, IERR )
*
* A = A[::-1]
*
      CALL SLAREV( 'R', N, N, A, LDA )
*
* A = A[p]
*
      CALL SLAPMR( .TRUE., N, N, A, LDA, JPVT )
*
* A = A * Q
*
      CALL SORMQR( 'R', 'N', N, N, N, B, LDB, WORK, A, LDA,
     $             WORK( N + 1 ), LWORK - N, IERR )
*
      IF( WANTQ ) THEN
*
* Initialize Q = J P J
* (J is the column index, I is the row index)
*
         DO J = 1, N
            DO I = 1, N
               Q( I, J ) = 0.0E0
            END DO
         END DO
         DO J = 1, N
            I = N - JPVT( N - J + 1 ) + 1
            Q( I, J ) = 1.0E0
         END DO
      END IF
*
      IF( WANTZ ) THEN
*
* Copy B into Z
*
         DO J = 1, N
            DO I = 1, N
               Z( I, J ) = B( I, J )
            END DO
         END DO
         CALL SORGQR( N, N, N, Z, LDZ, WORK, WORK( N + 1 ),
     $                LWORK - N, IERR )
         CALL SLAREV( 'C', N, N, Z, LDZ )
      END IF
*
* A = A[::-1][:, ::-1]
*
      CALL SLAREV( 'B', N, N, A, LDA )
*
* B = B.T[::-1][:, ::-1]
* (Zero strictly lower triangle to preserve the R matrix from SGEQP3.
* J is the column, I is the row.)
*
      DO J = 1, N - 1
         DO I = J + 1, N
            B( I, J ) = 0.0E0
         END DO
      END DO
*
      CALL SLATRN( N, B, LDB )
      CALL SLAREV( 'B', N, N, B, LDB )
*
* Store norm and tolerance in the first two elements of WORK
*
      WORK( 1 ) = NORMB
      WORK( 2 ) = TOL
*
      RETURN
*
* End of SGGPRP
*
      END