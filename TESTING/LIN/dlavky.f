*> \brief \b DLAVKY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAVKY( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B,
*                          LDB, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIAG, TRANS, UPLO
*       INTEGER            INFO, LDA, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAVKY  performs one of the matrix-vector operations
*>    x := A*x  or  x := A'*x,
*> where x is an N element vector and A is one of the factors
*> from the block U*D*U' or L*D*L' factorization computed by SKYTRF.
*>
*> If TRANS = 'N', multiplies by U  or U * D  (or L  or L * D)
*> If TRANS = 'T', multiplies by U' or D * U' (or L' or D * L')
*> If TRANS = 'C', multiplies by U' or D * U' (or L' or D * L')
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the factor stored in A is upper or lower
*>          triangular.
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>          Specifies the operation to be performed:
*>          = 'N':  x := A*x
*>          = 'T':  x := A'*x
*>          = 'C':  x := A'*x
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>          Specifies whether or not the diagonal blocks are unit
*>          matrices.  If the diagonal blocks are assumed to be unit,
*>          then A = U or A = L, otherwise A = U*D or A = L*D.
*>          = 'U':  Diagonal blocks are assumed to be unit matrices.
*>          = 'N':  Diagonal blocks are assumed to be non-unit matrices.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of rows and columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NRHS
*> \verbatim
*>          NRHS is INTEGER
*>          The number of right hand sides, i.e., the number of vectors
*>          x to be multiplied by A.  NRHS >= 0.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by SKYTRF.
*>          Stored as a 2-D triangular matrix.
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
*>          Details of the interchanges and the block structure of D,
*>          as determined by SKYTRF.
*>
*>          The elements of array IPIV are combined in pair, and the first 
*>          (if UPLO = 'U') or the second (if UPLO = 'L') element in
*>          the pair always keeps the value 0.  If N is odd, the first
*>          (if UPLO = 'U') or the last (if UPLO = 'L') element of IPIV is
*>          0, which is the only element not in pair. So we only use the
*>          first (if UPLO = 'L') or the second (if UPLO = 'U') element in
*>          the pair to determine the interchanges.
*>
*>          If IPIV(k)
*>          = 0: there was no interchange.
*>          > 0: rows and columns k-1 and IPIV(k) were interchanged, if
*>               UPLO = 'U', and rows and columns k+1 and IPIV(k) were
*>               interchanged, if UPLO = 'L'.
*>          < 0: rows and columns k and k-1 were interchanged,
*>               then rows and columns k-1 and -IPIV(k) were interchanged, if
*>               UPLO = 'U', and rows and columns k and k+1 were interchanged,
*>               then rows and columns k+1 and -IPIV(k) were interchanged, if
*>               UPLO = 'L'.
*> \endverbatim
*>
*> \param[in,out] B
*> \verbatim
*>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*>          On entry, B contains NRHS vectors of length N.
*>          On exit, B is overwritten with the product A * B.
*> \endverbatim
*>
*> \param[in] LDB
*> \verbatim
*>          LDB is INTEGER
*>          The leading dimension of the array B.  LDB >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -k, the k-th argument had an illegal value
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
*> \ingroup double_lin
*
*  =====================================================================
      SUBROUTINE DLAVKY( UPLO, TRANS, DIAG, N, NRHS, A, LDA, IPIV, B,
     $                   LDB, INFO )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, TRANS, UPLO
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT
      INTEGER            J, K, KP
      DOUBLE PRECISION   D11, D12, D21, D22, T1, T2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMV, DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.
     $         LSAME( TRANS, 'T' ) .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( DIAG, 'U' ) .AND. .NOT.LSAME( DIAG, 'N' ) )
     $          THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLAVKY ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*------------------------------------------
*
*     Compute  B := A * B  (No transpose)
*
*------------------------------------------
      IF( LSAME( TRANS, 'N' ) ) THEN
*
*        Compute  B := U*B
*        where U = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Loop forward applying the transformations.
*
            K = MOD(N, 2) + 1
   10       CONTINUE
            IF( K.GE.N )
     $         GO TO 30
*
*           2 x 2 pivot block
*
*           Multiply by the diagonal block if forming U * D.
*
            IF( NOUNIT ) THEN
               D11 = ZERO
               D22 = ZERO
               D12 = A( K, K+1 )
               D21 = -D12
               DO 20 J = 1, NRHS
                  T1 = B( K, J )
                  T2 = B( K+1, J )
                  B( K, J ) = D11*T1 + D12*T2
                  B( K+1, J ) = D21*T1 + D22*T2
   20          CONTINUE
            END IF
*
*           Multiply by  P(K) * inv(U(K))  if K > 1.
*
            IF( K.GT.1 ) THEN
*
*              Apply the transformations.
*
               CALL DGER( K-1, NRHS, ONE, A( 1, K ), 1, B( K, 1 ),
     $                    LDB, B( 1, 1 ), LDB )
               CALL DGER( K-1, NRHS, ONE, A( 1, K+1 ), 1,
     $                    B( K+1, 1 ), LDB, B( 1, 1 ), LDB )
*
*              Interchange if P(K) .ne. I.
*
               KP = IPIV( K+1 )
               IF( KP.GT.0 ) THEN
                  CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               ELSEIF( KP.LT.0 ) THEN
                  CALL DSWAP( NRHS, B( K, 1 ), LDB, B( -KP, 1 ), LDB )
                  CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K+1, 1 ), LDB )
               END IF
            END IF
            K = K + 2
            GO TO 10
   30       CONTINUE
*
*        Compute  B := L*B
*        where L = P(1)*inv(L(1))* ... *P(m)*inv(L(m)) .
*
         ELSE
*
*           Loop backward applying the transformations to B.
*
            K = N - MOD(N, 2)
   40       CONTINUE
            IF( K.LE.1 )
     $         GO TO 60
*
*           Test the pivot index.  A 2 x 2 pivot was used.
*
*           2 x 2 pivot block:
*
*           Multiply by the diagonal block if forming L * D.
*
            IF( NOUNIT ) THEN
               D11 = ZERO
               D22 = ZERO
               D21 = A( K, K-1 )
               D12 = -D21
               DO 50 J = 1, NRHS
                  T1 = B( K-1, J )
                  T2 = B( K, J )
                  B( K-1, J ) = D11*T1 + D12*T2
                  B( K, J ) = D21*T1 + D22*T2
   50          CONTINUE
            END IF
*
*           Multiply by  P(K) * inv(L(K))  if K < N.
*
            IF( K.NE.N ) THEN
*
*              Apply the transformation.
*
               CALL DGER( N-K, NRHS, ONE, A( K+1, K ), 1, B( K, 1 ),
     $                    LDB, B( K+1, 1 ), LDB )
               CALL DGER( N-K, NRHS, ONE, A( K+1, K-1 ), 1,
     $                    B( K-1, 1 ), LDB, B( K+1, 1 ), LDB )
*
*              Interchange if a permutation was applied at the
*              K-th step of the factorization.
*
               KP = IPIV( K-1 )
               IF( KP.GT.0 ) THEN
                  CALL DSWAP( NRHS, B( K, 1 ), LDB, B( KP, 1 ), LDB )
               ELSEIF( KP.LT.0 ) THEN
                  CALL DSWAP( NRHS, B( K, 1 ), LDB, B( -KP, 1 ), LDB )
                  CALL DSWAP( NRHS, B( K, 1 ), LDB, B( K-1, 1 ), LDB )
               END IF
            END IF
            K = K - 2
            GO TO 40
   60       CONTINUE
         END IF
*----------------------------------------
*
*     Compute  B := A' * B  (transpose)
*
*----------------------------------------
      ELSE
*
*        Form  B := U'*B
*        where U  = P(m)*inv(U(m))* ... *P(1)*inv(U(1))
*        and   U' = inv(U'(1))*P(1)* ... *inv(U'(m))*P(m)
*
         IF( LSAME( UPLO, 'U' ) ) THEN
*
*           Loop backward applying the transformations.
*
            K = N
   70       CONTINUE
            IF( K.LE.1 )
     $         GO TO 90
*
*           2 x 2 pivot block.
*
            IF( K.GT.2 ) THEN
*
*              Interchange if P(K) .ne. I.
*
               KP = IPIV( K )
               IF( KP.GT.0 ) THEN
                  CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( KP, 1 ),
     $                        LDB )
               ELSEIF( KP.LT.0 ) THEN
                  CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( K, 1 ),
     $                        LDB )
                  CALL DSWAP( NRHS, B( K-1, 1 ), LDB, B( -KP, 1 ),
     $                        LDB )
               ENDIF
*
*              Apply the transformations
*
               CALL DGEMV( 'Transpose', K-2, NRHS, ONE, B, LDB,
     $                     A( 1, K ), 1, ONE, B( K, 1 ), LDB )
               CALL DGEMV( 'Transpose', K-2, NRHS, ONE, B, LDB,
     $                     A( 1, K-1 ), 1, ONE, B( K-1, 1 ), LDB )
            END IF
*
*           Multiply by the diagonal block if non-unit.
*
            IF( NOUNIT ) THEN
               D11 = ZERO
               D22 = ZERO
               D12 = A( K-1, K )
               D21 = -D12
               DO 80 J = 1, NRHS
                  T1 = B( K-1, J )
                  T2 = B( K, J )
                  B( K-1, J ) = D11*T1 + D12*T2
                  B( K, J ) = D21*T1 + D22*T2
   80          CONTINUE
            END IF
            K = K - 2
            GO TO 70
   90       CONTINUE
*
*        Form  B := L'*B
*        where L  = P(1)*inv(L(1))* ... *P(m)*inv(L(m))
*        and   L' = inv(L'(m))*P(m)* ... *inv(L'(1))*P(1)
*
         ELSE
*
*           Loop forward applying the L-transformations.
*
            K = 1
  100       CONTINUE
            IF( K.GE.N )
     $         GO TO 120
*
*           2 x 2 pivot block
*
            IF( K.LT.N-1 ) THEN
*
*           Interchange if P(K) .ne. I.
*
               KP = IPIV( K )
               IF( KP.GT.0 ) THEN
                  CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( KP, 1 ),
     $                        LDB )
               ELSEIF( KP.LT.0 ) THEN
                  CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( K, 1 ),
     $                        LDB )
                  CALL DSWAP( NRHS, B( K+1, 1 ), LDB, B( -KP, 1 ),
     $                        LDB )
               ENDIF
*
*              Apply the transformation
*
               CALL DGEMV( 'Transpose', N-K-1, NRHS, ONE,
     $                     B( K+2, 1 ), LDB, A( K+2, K+1 ), 1, ONE,
     $                     B( K+1, 1 ), LDB )
               CALL DGEMV( 'Transpose', N-K-1, NRHS, ONE,
     $                     B( K+2, 1 ), LDB, A( K+2, K ), 1, ONE,
     $                     B( K, 1 ), LDB )
            END IF
*
*           Multiply by the diagonal block if non-unit.
*
            IF( NOUNIT ) THEN
               D11 = ZERO
               D22 = ZERO
               D21 = A( K+1, K )
               D12 = -D21
               DO 110 J = 1, NRHS
                  T1 = B( K, J )
                  T2 = B( K+1, J )
                  B( K, J ) = D11*T1 + D12*T2
                  B( K+1, J ) = D21*T1 + D22*T2
  110          CONTINUE
            END IF
            K = K + 2
            GO TO 100
  120       CONTINUE
         END IF
*
      END IF
      RETURN
*
*     End of DLAVKY
*
      END
