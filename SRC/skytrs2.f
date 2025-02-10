*> \brief \b SKYTRS2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SKYTRS2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skytrs2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skytrs2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skytrs2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
*                           WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKYTRS2 solves a system of linear equations A*X = B with a real
*> skew-symmetric matrix A using the factorization A = U*D*U**T or
*> A = L*D*L**T computed by SKYTRF and converted by SKYCONV.
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
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
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
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by SKYTRF.
*>          Note that A is input / output. This might be counter-intuitive,
*>          and one may think that A is input only. A is input / output. This
*>          is because, at the start of the subroutine, we permute A in a
*>          "better" form and then we permute A back to its original form at
*>          the end.
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
*>          Details of the interchanges and the block structure of D
*>          as determined by SKYTRF.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is REAL array, dimension (LDB,NRHS)
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
*>          WORK is REAL array, dimension (N)
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
*> \ingroup kytrs2
*
*  =====================================================================
      SUBROUTINE SKYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
     $                    WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, IINFO, K, KP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           SSCAL, SKYCONV, SSWAP, STRSM,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
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
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKYTRS2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
*     Convert A
*
      CALL SKYCONV( UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO )
*
      IF( UPPER ) THEN
*
*        Solve A*X = B, where A = U*D*U**T.
*
*       P**T * B
        K=N
        DO WHILE ( K .GE. 2 )
         IF( IPIV( K ).GT.0 ) THEN
*           2 x 2 diagonal block
*           Interchange rows K-1 and IPIV(K).
            KP = IPIV( K )
            CALL SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K), then K and K-1.
            KP = -IPIV( K )
            CALL SSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
            CALL SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
         END IF
         K=K-2
        END DO
*
*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
*
        CALL STRSM('L','U','N','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*  Compute D \ B -> B   [ D \ (U \P**T * B) ]
*
         I=N
         DO WHILE ( I .GE. 2 )
            CALL SSCAL( NRHS, -ONE / WORK( I ), B( I, 1 ), LDB )
            CALL SSCAL( NRHS, ONE / WORK( I ), B( I-1, 1 ), LDB )
            CALL SSWAP( NRHS, B( I, 1 ), LDB, B( I-1, 1 ), LDB )
            I = I - 2
         END DO
*
*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
*
         CALL STRSM('L','U','T','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
*
        K=2
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
*           2 x 2 diagonal block
*           Interchange rows K-1 and IPIV(K).
            KP = IPIV( K )
            CALL SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )  
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K and K-1, then K-1 and -IPIV(K).
            KP = -IPIV( K )
            CALL SSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            CALL SSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
         ENDIF
         K=K+2
        END DO
*
      ELSE
*
*        Solve A*X = B, where A = L*D*L**T.
*
*       P**T * B
        K=1
        DO WHILE ( K .LE. N-1 )
         IF( IPIV( K ).GT.0 ) THEN
*           2 x 2 diagonal block
*           Interchange rows K+1 and IPIV(K).
            KP = IPIV( K )
            CALL SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K+1 and -IPIV(K), then K and K+1.
            KP = -IPIV( K )
            CALL SSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
            CALL SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
         ENDIF
         K=K+2
        END DO
*
*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
*
        CALL STRSM('L','L','N','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*  Compute D \ B -> B   [ D \ (L \P**T * B) ]
*
         I=1
         DO WHILE ( I .LE. N-1 )
            CALL SSCAL( NRHS, -ONE / WORK( I ), B( I, 1 ), LDB )
            CALL SSCAL( NRHS, ONE / WORK( I ), B( I+1, 1 ), LDB )
            CALL SSWAP( NRHS, B( I, 1 ), LDB, B( I+1, 1 ), LDB )
            I = I + 2
         END DO
*
*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
*
        CALL STRSM('L','L','T','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
*
        K=N-1
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
*           2 x 2 diagonal block
*           Interchange rows K+1 and IPIV(K).
            KP = IPIV( K )
            CALL SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K and K+1, then K+1 and -IPIV(K).
            KP = -IPIV( K )
            CALL SSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            CALL SSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
         ENDIF
         K=K-2
        END DO
*
      END IF
*
*     Revert A
*
      CALL SKYCONV( UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO )
*
      RETURN
*
*     End of SKYTRS2
*
      END
