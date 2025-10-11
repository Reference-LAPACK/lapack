*> \brief \b DKYTRS2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DKYTRS2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dkytrs2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dkytrs2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dkytrs2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DKYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
*                           WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DKYTRS2 solves a system of linear equations A*X = B with a real
*> skew-symmetric matrix A using the factorization A = U*D*U**T or
*> A = L*D*L**T computed by SKYTRF and converted by DKYCONV.
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
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
*>          WORK is DOUBLE PRECISION array, dimension (N)
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
*> \ingroup hetrs2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Derived from subroutine dsytrs2.
*>
*> \endverbatim
*
*> \par Contributors:
*  ==================
*>
*>    Shuo Zheng, China, Jul 2025 \n
*>
*  =====================================================================
      SUBROUTINE DKYTRS2( UPLO, N, NRHS, A, LDA, IPIV, B, LDB,
     $                    WORK, INFO )
      IMPLICIT NONE
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
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
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
      EXTERNAL           DSCAL, DKYCONV, DSWAP, DTRSM,
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
         CALL XERBLA( 'DKYTRS2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 .OR. (MOD(N,2).NE.0) )
     $   RETURN
*
*     Convert A
*
      CALL DKYCONV( UPLO, 'C', N, A, LDA, IPIV, WORK, IINFO )
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
            CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K-1 and -IPIV(K), then K and K-1.
            KP = -IPIV( K )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
         END IF
         K=K-2
        END DO
*
*  Compute (U \P**T * B) -> B    [ (U \P**T * B) ]
*
        CALL DTRSM('L','U','N','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*  Compute D \ B -> B   [ D \ (U \P**T * B) ]
*
         I=N
         DO WHILE ( I .GE. 2 )
            CALL DSCAL( NRHS, -ONE / WORK( I ), B( I, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / WORK( I ), B( I-1, 1 ), LDB )
            CALL DSWAP( NRHS, B( I, 1 ), LDB, B( I-1, 1 ), LDB )
            I = I - 2
         END DO
*
*      Compute (U**T \ B) -> B   [ U**T \ (D \ (U \P**T * B) ) ]
*
         CALL DTRSM('L','U','T','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*       P * B  [ P * (U**T \ (D \ (U \P**T * B) )) ]
*
        K=2
        DO WHILE ( K .LE. N )
         IF( IPIV( K ).GT.0 ) THEN
*           2 x 2 diagonal block
*           Interchange rows K-1 and IPIV(K).
            KP = IPIV( K )
            CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )  
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K and K-1, then K-1 and -IPIV(K).
            KP = -IPIV( K )
            CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
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
            CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K+1 and -IPIV(K), then K and K+1.
            KP = -IPIV( K )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
            CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
         ENDIF
         K=K+2
        END DO
*
*  Compute (L \P**T * B) -> B    [ (L \P**T * B) ]
*
        CALL DTRSM('L','L','N','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*  Compute D \ B -> B   [ D \ (L \P**T * B) ]
*
         I=1
         DO WHILE ( I .LE. N-1 )
            CALL DSCAL( NRHS, -ONE / WORK( I ), B( I, 1 ), LDB )
            CALL DSCAL( NRHS, ONE / WORK( I ), B( I+1, 1 ), LDB )
            CALL DSWAP( NRHS, B( I, 1 ), LDB, B( I+1, 1 ), LDB )
            I = I + 2
         END DO
*
*  Compute (L**T \ B) -> B   [ L**T \ (D \ (L \P**T * B) ) ]
*
        CALL DTRSM('L','L','T','U',N,NRHS,ONE,A,LDA,B,LDB)
*
*       P * B  [ P * (L**T \ (D \ (L \P**T * B) )) ]
*
        K=N-1
        DO WHILE ( K .GE. 1 )
         IF( IPIV( K ).GT.0 ) THEN
*           2 x 2 diagonal block
*           Interchange rows K+1 and IPIV(K).
            KP = IPIV( K )
            CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
         ELSEIF ( IPIV( K ).LT.0) THEN
*           2 x 2 diagonal block
*           Interchange rows K and K+1, then K+1 and -IPIV(K).
            KP = -IPIV( K )
            CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ), LDB )
            CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
         ENDIF
         K=K-2
        END DO
*
      END IF
*
*     Revert A
*
      CALL DKYCONV( UPLO, 'R', N, A, LDA, IPIV, WORK, IINFO )
*
      RETURN
*
*     End of DKYTRS2
*
      END
