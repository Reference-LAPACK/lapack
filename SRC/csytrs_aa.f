*> \brief \b CSYTRS_AA
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CSYTRS_AA + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/csytrs_aa.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/csytrs_aa.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/csytrs_aa.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
*                             WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CSYTRS_AA solves a system of linear equations A*X = B with a complex
*> symmetric matrix A using the factorization A = U**T*T*U or
*> A = L*T*L**T computed by CSYTRF_AA.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the details of the factorization are stored
*>          as an upper or lower triangular matrix.
*>          = 'U':  Upper triangular, form is A = U**T*T*U;
*>          = 'L':  Lower triangular, form is A = L*T*L**T.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of columns
*>          of the matrix B.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is COMPLEX array, dimension (LDA,N)
*>          Details of factors computed by CSYTRF_AA.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges as computed by CSYTRF_AA.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is COMPLEX array, dimension (LDB,NRHS)
*>          On entry, the right hand side matrix B.
*>          On exit, the solution matrix X.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,3*N-2).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
*> \ingroup hetrs_aa
*
*  =====================================================================
      SUBROUTINE CSYTRS_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
     $                      WORK, LWORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
      IMPLICIT NONE
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            N, NRHS, LDA, LDB, LWORK, INFO
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
      COMPLEX            ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            K, KP, LWKOPT
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SROUNDUP_LWORK
      EXTERNAL           LSAME, SROUNDUP_LWORK
*     ..
*     .. External Subroutines ..
      EXTERNAL           CLACPY, CGTSV, CSWAP, CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, 3*N-2 ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CSYTRS_AA', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         LWKOPT = (3*N-2)
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U**T*T*U.
*
*        1) Forward substitution with U**T
*
         IF( N.GT.1 ) THEN
*
*           Pivot, P**T * B -> B
*
            DO K = 1, N
               KP = IPIV( K )
               IF( KP.NE.K )
     $            CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO
*
*           Compute U**T \ B -> B    [ (U**T \P**T * B) ]
*
            CALL CTRSM( 'L', 'U', 'T', 'U', N-1, NRHS, ONE, A( 1,
     $                  2 ),
     $                  LDA, B( 2, 1 ), LDB)
         END IF
*
*        2) Solve with triangular matrix T
*
*        Compute T \ B -> B   [ T \ (U**T \P**T * B) ]
*
         CALL CLACPY( 'F', 1, N, A( 1, 1 ), LDA+1, WORK( N ), 1)
         IF( N.GT.1 ) THEN
            CALL CLACPY( 'F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 1 ),
     $                   1 )
            CALL CLACPY( 'F', 1, N-1, A( 1, 2 ), LDA+1, WORK( 2*N ),
     $                   1 )
         END IF
         CALL CGTSV( N, NRHS, WORK( 1 ), WORK( N ), WORK( 2*N ), B,
     $               LDB,
     $               INFO )
*
*        3) Backward substitution with U
*
         IF( N.GT.1 ) THEN
*
*           Compute U \ B -> B   [ U \ (T \ (U**T \P**T * B) ) ]
*
            CALL CTRSM( 'L', 'U', 'N', 'U', N-1, NRHS, ONE, A( 1,
     $                  2 ),
     $                  LDA, B( 2, 1 ), LDB)
*
*           Pivot, P * B -> B  [ P * (U**T \ (T \ (U \P**T * B) )) ]
*
            DO K = N, 1, -1
               KP = IPIV( K )
               IF( KP.NE.K )
     $            CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO
         END IF
*
      ELSE
*
*        Solve A*X = B, where A = L*T*L**T.
*
*        1) Forward substitution with L
*
         IF( N.GT.1 ) THEN
*
*           Pivot, P**T * B -> B
*
            DO K = 1, N
               KP = IPIV( K )
               IF( KP.NE.K )
     $            CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO
*
*           Compute L \ B -> B    [ (L \P**T * B) ]
*
            CALL CTRSM( 'L', 'L', 'N', 'U', N-1, NRHS, ONE, A( 2,
     $                  1 ),
     $                  LDA, B( 2, 1 ), LDB)
         END IF
*
*        2) Solve with triangular matrix T
*
*
*        Compute T \ B -> B   [ T \ (L \P**T * B) ]
*
         CALL CLACPY( 'F', 1, N, A(1, 1), LDA+1, WORK(N), 1)
         IF( N.GT.1 ) THEN
            CALL CLACPY( 'F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 1 ),
     $                   1 )
            CALL CLACPY( 'F', 1, N-1, A( 2, 1 ), LDA+1, WORK( 2*N ),
     $                   1 )
         END IF
         CALL CGTSV( N, NRHS, WORK( 1 ), WORK(N), WORK( 2*N ), B,
     $               LDB,
     $               INFO)
*
*        3) Backward substitution with L**T
*
         IF( N.GT.1 ) THEN
*
*           Compute (L**T \ B) -> B   [ L**T \ (T \ (L \P**T * B) ) ]
*
            CALL CTRSM( 'L', 'L', 'T', 'U', N-1, NRHS, ONE, A( 2,
     $                  1 ),
     $                  LDA, B( 2, 1 ), LDB)
*
*           Pivot, P * B -> B  [ P * (L**T \ (T \ (L \P**T * B) )) ]
*
            DO K = N, 1, -1
               KP = IPIV( K )
               IF( KP.NE.K )
     $            CALL CSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
            END DO
         END IF
*
      END IF
*
      RETURN
*
*     End of CSYTRS_AA
*
      END
