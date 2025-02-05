*> \brief \b ZHPGST
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZHPGST + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhpgst.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhpgst.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhpgst.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, ITYPE, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         AP( * ), BP( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZHPGST reduces a complex Hermitian-definite generalized
*> eigenproblem to standard form, using packed storage.
*>
*> If ITYPE = 1, the problem is A*x = lambda*B*x,
*> and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
*>
*> If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
*> B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
*>
*> B must have been previously factorized as U**H*U or L*L**H by ZPPTRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ITYPE
*> \verbatim
*>          ITYPE is INTEGER
*>          = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
*>          = 2 or 3: compute U*A*U**H or L**H*A*L.
*> \endverbatim
*>
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          = 'U':  Upper triangle of A is stored and B is factored as
*>                  U**H*U;
*>          = 'L':  Lower triangle of A is stored and B is factored as
*>                  L*L**H.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrices A and B.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] AP
*> \verbatim
*>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
*>          On entry, the upper or lower triangle of the Hermitian matrix
*>          A, packed columnwise in a linear array.  The j-th column of A
*>          is stored in the array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
*>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
*>
*>          On exit, if INFO = 0, the transformed matrix, stored in the
*>          same format as A.
*> \endverbatim
*>
*> \param[in] BP
*> \verbatim
*>          BP is COMPLEX*16 array, dimension (N*(N+1)/2)
*>          The triangular factor from the Cholesky factorization of B,
*>          stored in the same format as A, as returned by ZPPTRF.
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
*> \ingroup hpgst
*
*  =====================================================================
      SUBROUTINE ZHPGST( ITYPE, UPLO, N, AP, BP, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, ITYPE, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         AP( * ), BP( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, HALF
      PARAMETER          ( ONE = 1.0D+0, HALF = 0.5D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, J1, J1J1, JJ, K, K1, K1K1, KK
      DOUBLE PRECISION   AJJ, AKK, BJJ, BKK
      COMPLEX*16         CT
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZDSCAL, ZHPMV, ZHPR2,
     $                   ZTPMV,
     $                   ZTPSV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( ITYPE.LT.1 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHPGST', -INFO )
         RETURN
      END IF
*
      IF( ITYPE.EQ.1 ) THEN
         IF( UPPER ) THEN
*
*           Compute inv(U**H)*A*inv(U)
*
*           J1 and JJ are the indices of A(1,j) and A(j,j)
*
            JJ = 0
            DO 10 J = 1, N
               J1 = JJ + 1
               JJ = JJ + J
*
*              Compute the j-th column of the upper triangle of A
*
               AP( JJ ) = DBLE( AP( JJ ) )
               BJJ = DBLE( BP( JJ ) )
               CALL ZTPSV( UPLO, 'Conjugate transpose', 'Non-unit',
     $                     J,
     $                     BP, AP( J1 ), 1 )
               CALL ZHPMV( UPLO, J-1, -CONE, AP, BP( J1 ), 1, CONE,
     $                     AP( J1 ), 1 )
               CALL ZDSCAL( J-1, ONE / BJJ, AP( J1 ), 1 )
               AP( JJ ) = ( AP( JJ )-ZDOTC( J-1, AP( J1 ), 1,
     $             BP( J1 ),
     $                    1 ) ) / BJJ
   10       CONTINUE
         ELSE
*
*           Compute inv(L)*A*inv(L**H)
*
*           KK and K1K1 are the indices of A(k,k) and A(k+1,k+1)
*
            KK = 1
            DO 20 K = 1, N
               K1K1 = KK + N - K + 1
*
*              Update the lower triangle of A(k:n,k:n)
*
               AKK = DBLE( AP( KK ) )
               BKK = DBLE( BP( KK ) )
               AKK = AKK / BKK**2
               AP( KK ) = AKK
               IF( K.LT.N ) THEN
                  CALL ZDSCAL( N-K, ONE / BKK, AP( KK+1 ), 1 )
                  CT = -HALF*AKK
                  CALL ZAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL ZHPR2( UPLO, N-K, -CONE, AP( KK+1 ), 1,
     $                        BP( KK+1 ), 1, AP( K1K1 ) )
                  CALL ZAXPY( N-K, CT, BP( KK+1 ), 1, AP( KK+1 ), 1 )
                  CALL ZTPSV( UPLO, 'No transpose', 'Non-unit', N-K,
     $                        BP( K1K1 ), AP( KK+1 ), 1 )
               END IF
               KK = K1K1
   20       CONTINUE
         END IF
      ELSE
         IF( UPPER ) THEN
*
*           Compute U*A*U**H
*
*           K1 and KK are the indices of A(1,k) and A(k,k)
*
            KK = 0
            DO 30 K = 1, N
               K1 = KK + 1
               KK = KK + K
*
*              Update the upper triangle of A(1:k,1:k)
*
               AKK = DBLE( AP( KK ) )
               BKK = DBLE( BP( KK ) )
               CALL ZTPMV( UPLO, 'No transpose', 'Non-unit', K-1, BP,
     $                     AP( K1 ), 1 )
               CT = HALF*AKK
               CALL ZAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL ZHPR2( UPLO, K-1, CONE, AP( K1 ), 1, BP( K1 ), 1,
     $                     AP )
               CALL ZAXPY( K-1, CT, BP( K1 ), 1, AP( K1 ), 1 )
               CALL ZDSCAL( K-1, BKK, AP( K1 ), 1 )
               AP( KK ) = AKK*BKK**2
   30       CONTINUE
         ELSE
*
*           Compute L**H *A*L
*
*           JJ and J1J1 are the indices of A(j,j) and A(j+1,j+1)
*
            JJ = 1
            DO 40 J = 1, N
               J1J1 = JJ + N - J + 1
*
*              Compute the j-th column of the lower triangle of A
*
               AJJ = DBLE( AP( JJ ) )
               BJJ = DBLE( BP( JJ ) )
               AP( JJ ) = AJJ*BJJ + ZDOTC( N-J, AP( JJ+1 ), 1,
     $                    BP( JJ+1 ), 1 )
               CALL ZDSCAL( N-J, BJJ, AP( JJ+1 ), 1 )
               CALL ZHPMV( UPLO, N-J, CONE, AP( J1J1 ), BP( JJ+1 ),
     $                     1,
     $                     CONE, AP( JJ+1 ), 1 )
               CALL ZTPMV( UPLO, 'Conjugate transpose', 'Non-unit',
     $                     N-J+1, BP( JJ ), AP( JJ ), 1 )
               JJ = J1J1
   40       CONTINUE
         END IF
      END IF
      RETURN
*
*     End of ZHPGST
*
      END
