*> \brief \b CHPTRI
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CHPTRI + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chptri.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chptri.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chptri.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       COMPLEX            AP( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CHPTRI computes the inverse of a complex Hermitian indefinite matrix
*> A in packed storage using the factorization A = U*D*U**H or
*> A = L*D*L**H computed by CHPTRF.
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
*>          = 'U':  Upper triangular, form is A = U*D*U**H;
*>          = 'L':  Lower triangular, form is A = L*D*L**H.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] AP
*> \verbatim
*>          AP is COMPLEX array, dimension (N*(N+1)/2)
*>          On entry, the block diagonal matrix D and the multipliers
*>          used to obtain the factor U or L as computed by CHPTRF,
*>          stored as a packed triangular matrix.
*>
*>          On exit, if INFO = 0, the (Hermitian) inverse of the original
*>          matrix, stored as a packed triangular matrix. The j-th column
*>          of inv(A) is stored in the array AP as follows:
*>          if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;
*>          if UPLO = 'L',
*>             AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by CHPTRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (N)
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its
*>               inverse could not be computed.
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
*> \ingroup hptri
*
*  =====================================================================
      SUBROUTINE CHPTRI( UPLO, N, AP, IPIV, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX            AP( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      COMPLEX            CONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ),
     $                   ZERO = ( 0.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, K, KC, KCNEXT, KP, KPC, KSTEP, KX, NPP
      REAL               AK, AKP1, D, T
      COMPLEX            AKKP1, TEMP
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX            CDOTC
      EXTERNAL           LSAME, CDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           CCOPY, CHPMV, CSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, REAL
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CHPTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check that the diagonal matrix D is nonsingular.
*
      IF( UPPER ) THEN
*
*        Upper triangular storage: examine D from bottom to top
*
         KP = N*( N+1 ) / 2
         DO 10 INFO = N, 1, -1
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO )
     $         RETURN
            KP = KP - INFO
   10    CONTINUE
      ELSE
*
*        Lower triangular storage: examine D from top to bottom.
*
         KP = 1
         DO 20 INFO = 1, N
            IF( IPIV( INFO ).GT.0 .AND. AP( KP ).EQ.ZERO )
     $         RETURN
            KP = KP + N - INFO + 1
   20    CONTINUE
      END IF
      INFO = 0
*
      IF( UPPER ) THEN
*
*        Compute inv(A) from the factorization A = U*D*U**H.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         K = 1
         KC = 1
   30    CONTINUE
*
*        If K > N, exit from loop.
*
         IF( K.GT.N )
     $      GO TO 50
*
         KCNEXT = KC + K
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            AP( KC+K-1 ) = ONE / REAL( AP( KC+K-1 ) )
*
*           Compute column K of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL CCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL CHPMV( UPLO, K-1, -CONE, AP, WORK, 1, ZERO,
     $                     AP( KC ), 1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) -
     $                        REAL( CDOTC( K-1, WORK, 1, AP( KC ),
     $                              1 ) )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( AP( KCNEXT+K-1 ) )
            AK = REAL( AP( KC+K-1 ) ) / T
            AKP1 = REAL( AP( KCNEXT+K ) ) / T
            AKKP1 = AP( KCNEXT+K-1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KC+K-1 ) = AKP1 / D
            AP( KCNEXT+K ) = AK / D
            AP( KCNEXT+K-1 ) = -AKKP1 / D
*
*           Compute columns K and K+1 of the inverse.
*
            IF( K.GT.1 ) THEN
               CALL CCOPY( K-1, AP( KC ), 1, WORK, 1 )
               CALL CHPMV( UPLO, K-1, -CONE, AP, WORK, 1, ZERO,
     $                     AP( KC ), 1 )
               AP( KC+K-1 ) = AP( KC+K-1 ) -
     $                        REAL( CDOTC( K-1, WORK, 1, AP( KC ),
     $                              1 ) )
               AP( KCNEXT+K-1 ) = AP( KCNEXT+K-1 ) -
     $                            CDOTC( K-1, AP( KC ), 1,
     $                                   AP( KCNEXT ),
     $                            1 )
               CALL CCOPY( K-1, AP( KCNEXT ), 1, WORK, 1 )
               CALL CHPMV( UPLO, K-1, -CONE, AP, WORK, 1, ZERO,
     $                     AP( KCNEXT ), 1 )
               AP( KCNEXT+K ) = AP( KCNEXT+K ) -
     $                          REAL( CDOTC( K-1, WORK, 1,
     $                                AP( KCNEXT ),
     $                          1 ) )
            END IF
            KSTEP = 2
            KCNEXT = KCNEXT + K + 1
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the leading
*           submatrix A(1:k+1,1:k+1)
*
            KPC = ( KP-1 )*KP / 2 + 1
            CALL CSWAP( KP-1, AP( KC ), 1, AP( KPC ), 1 )
            KX = KPC + KP - 1
            DO 40 J = KP + 1, K - 1
               KX = KX + J - 1
               TEMP = CONJG( AP( KC+J-1 ) )
               AP( KC+J-1 ) = CONJG( AP( KX ) )
               AP( KX ) = TEMP
   40       CONTINUE
            AP( KC+KP-1 ) = CONJG( AP( KC+KP-1 ) )
            TEMP = AP( KC+K-1 )
            AP( KC+K-1 ) = AP( KPC+KP-1 )
            AP( KPC+KP-1 ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = AP( KC+K+K-1 )
               AP( KC+K+K-1 ) = AP( KC+K+KP-1 )
               AP( KC+K+KP-1 ) = TEMP
            END IF
         END IF
*
         K = K + KSTEP
         KC = KCNEXT
         GO TO 30
   50    CONTINUE
*
      ELSE
*
*        Compute inv(A) from the factorization A = L*D*L**H.
*
*        K is the main loop index, increasing from 1 to N in steps of
*        1 or 2, depending on the size of the diagonal blocks.
*
         NPP = N*( N+1 ) / 2
         K = N
         KC = NPP
   60    CONTINUE
*
*        If K < 1, exit from loop.
*
         IF( K.LT.1 )
     $      GO TO 80
*
         KCNEXT = KC - ( N-K+2 )
         IF( IPIV( K ).GT.0 ) THEN
*
*           1 x 1 diagonal block
*
*           Invert the diagonal block.
*
            AP( KC ) = ONE / REAL( AP( KC ) )
*
*           Compute column K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL CCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL CHPMV( UPLO, N-K, -CONE, AP( KC+N-K+1 ), WORK, 1,
     $                     ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - REAL( CDOTC( N-K, WORK, 1,
     $                    AP( KC+1 ), 1 ) )
            END IF
            KSTEP = 1
         ELSE
*
*           2 x 2 diagonal block
*
*           Invert the diagonal block.
*
            T = ABS( AP( KCNEXT+1 ) )
            AK = REAL( AP( KCNEXT ) ) / T
            AKP1 = REAL( AP( KC ) ) / T
            AKKP1 = AP( KCNEXT+1 ) / T
            D = T*( AK*AKP1-ONE )
            AP( KCNEXT ) = AKP1 / D
            AP( KC ) = AK / D
            AP( KCNEXT+1 ) = -AKKP1 / D
*
*           Compute columns K-1 and K of the inverse.
*
            IF( K.LT.N ) THEN
               CALL CCOPY( N-K, AP( KC+1 ), 1, WORK, 1 )
               CALL CHPMV( UPLO, N-K, -CONE, AP( KC+( N-K+1 ) ),
     $                     WORK,
     $                     1, ZERO, AP( KC+1 ), 1 )
               AP( KC ) = AP( KC ) - REAL( CDOTC( N-K, WORK, 1,
     $                    AP( KC+1 ), 1 ) )
               AP( KCNEXT+1 ) = AP( KCNEXT+1 ) -
     $                          CDOTC( N-K, AP( KC+1 ), 1,
     $                          AP( KCNEXT+2 ), 1 )
               CALL CCOPY( N-K, AP( KCNEXT+2 ), 1, WORK, 1 )
               CALL CHPMV( UPLO, N-K, -CONE, AP( KC+( N-K+1 ) ),
     $                     WORK,
     $                     1, ZERO, AP( KCNEXT+2 ), 1 )
               AP( KCNEXT ) = AP( KCNEXT ) -
     $                        REAL( CDOTC( N-K, WORK, 1,
     $                              AP( KCNEXT+2 ),
     $                        1 ) )
            END IF
            KSTEP = 2
            KCNEXT = KCNEXT - ( N-K+3 )
         END IF
*
         KP = ABS( IPIV( K ) )
         IF( KP.NE.K ) THEN
*
*           Interchange rows and columns K and KP in the trailing
*           submatrix A(k-1:n,k-1:n)
*
            KPC = NPP - ( N-KP+1 )*( N-KP+2 ) / 2 + 1
            IF( KP.LT.N )
     $         CALL CSWAP( N-KP, AP( KC+KP-K+1 ), 1, AP( KPC+1 ), 1 )
            KX = KC + KP - K
            DO 70 J = K + 1, KP - 1
               KX = KX + N - J + 1
               TEMP = CONJG( AP( KC+J-K ) )
               AP( KC+J-K ) = CONJG( AP( KX ) )
               AP( KX ) = TEMP
   70       CONTINUE
            AP( KC+KP-K ) = CONJG( AP( KC+KP-K ) )
            TEMP = AP( KC )
            AP( KC ) = AP( KPC )
            AP( KPC ) = TEMP
            IF( KSTEP.EQ.2 ) THEN
               TEMP = AP( KC-N+K-1 )
               AP( KC-N+K-1 ) = AP( KC-N+KP-1 )
               AP( KC-N+KP-1 ) = TEMP
            END IF
         END IF
*
         K = K - KSTEP
         KC = KCNEXT
         GO TO 60
   80    CONTINUE
      END IF
*
      RETURN
*
*     End of CHPTRI
*
      END
