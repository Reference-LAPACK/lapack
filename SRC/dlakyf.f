*> \brief \b DLAKYF computes a partial factorization of a real skew-symmetric matrix using the Bunch partial pivoting method.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download SLASYF + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlakyf.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlakyf.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlakyf.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLAKYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO
*       INTEGER            INFO, KB, LDA, LDW, N, NB
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLAKYF computes a partial factorization of a real skew-symmetric matrix A
*> using the Bunch partial pivoting method. The partial factorization has
*> the form:
*>
*> A  =  ( I  U12 ) ( A11  0  ) (  I       0       )  if UPLO = 'U', or:
*>       ( 0  U22 ) (  0   D  ) (  U12**T  U22**T )
*>
*> A  =  ( L11  0 ) (  D   0  ) (  L11**T  L21**T )  if UPLO = 'L'
*>       ( L21  I ) (  0  A22 ) (  0       I       )
*>
*> where the order of D is at most NB. The actual order is returned in the
*> argument KB, and is either NB or NB-1, or N if N <= NB.
*>
*> DLAKYF is an auxiliary routine called by DKYTRF. It uses blocked code
*> (calling Level 3 BLAS) to update the submatrix A11 (if UPLO = 'U') or
*> A22 (if UPLO = 'L').
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>          Specifies whether the upper or lower triangular part of the
*>          skew-symmetric matrix A is stored:
*>          = 'U':  Upper triangular
*>          = 'L':  Lower triangular
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The maximum number of columns of the matrix A that should be
*>          factored.  NB should be at least 2 to allow for 2-by-2 pivot
*>          blocks.
*> \endverbatim
*>
*> \param[out] KB
*> \verbatim
*>          KB is INTEGER
*>          The number of columns of A that were actually factored.
*>          KB is either NB-1 or NB, or N if N <= NB.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the skew-symmetric matrix A.  If UPLO = 'U', the
*>          strictly upper triangular part of A contains the upper
*>          triangular part of the matrix A, and the leading N-by-N lower
*>          triangular part of A is not referenced.  If UPLO = 'L', the
*>          strictly lower triangular part of A contains the lower
*>          triangular part of the matrix A, and the leading N-by-N upper
*>          triangular part of A is not referenced.
*>          On exit, A contains details of the partial factorization.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[out] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D.
*>
*>          If UPLO = 'U':
*>             Only the last KB elements of IPIV are set.
*>
*>             The elements of array IPIV are combined in pair, and the first 
*>             element in the pair always keeps the value 0.  If N is odd, the
*>             first element of IPIV is 0, which is the only element not in pair.
*>             So we only use the second element in the pair to determine the
*>             interchanges.
*>
*>             If IPIV(k)
*>             = 0: there was no interchange.
*>             > 0: rows and columns k-1 and IPIV(k) were interchanged.
*>             < 0: rows and columns k and k-1 were interchanged,
*>                  then rows and columns k-1 and -IPIV(k) were interchanged.
*>
*>          If UPLO = 'L':
*>             Only the first KB elements of IPIV are set.
*>
*>             The elements of array IPIV are combined in pair, and the second
*>             element in the pair always keeps the value 0.  If N is odd, the
*>             last element of IPIV is 0, which is the only element not in pair.
*>             So we only use the first element in the pair to determine the
*>             interchanges.
*>
*>             If IPIV(k)
*>             = 0: there was no interchange.
*>             > 0: rows and columns k+1 and IPIV(k) were interchanged。
*>             < 0: rows and columns k and k+1 were interchanged,
*>                  then rows and columns k+1 and -IPIV(k) were interchanged.
*> \endverbatim
*>
*> \param[out] W
*> \verbatim
*>          W is DOUBLE PRECISION array, dimension (LDW,NB)
*> \endverbatim
*>
*> \param[in] LDW
*> \verbatim
*>          LDW is INTEGER
*>          The leading dimension of the array W.  LDW >= max(1,N).
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
*>          > 0: if INFO = i, D(i-1,i) (if UPLO = 'U') or D(i+1,i) (if UPLO = 'L')
*>               is exactly zero.  The factorization has been completed,
*>               but the block diagonal matrix D is exactly singular,
*>               so the solution could not be computed.
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
*> \ingroup lahef
*
*> \par Contributors:
*  ==================
*>
*> \verbatim
*>
*>    Shuo Zheng, China, Jul 2025 \n
*>
*> \endverbatim
*
*  =====================================================================
      SUBROUTINE DLAKYF( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW,
     $					 INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IMAX1, IMAX2, J, JB, JJ, JMAX, JP, K,
     $                   KP, KW, KADJ
      DOUBLE PRECISION   ABSAKP1K, COLMAX1, COLMAX2
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IDAMAX
      EXTERNAL           LSAME, IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      KADJ = 0

*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Factorize the leading columns of A using the upper triangle
*        of A and working forwards, and compute the matrix W = U12*D
*        for use in updating A11
*
*        K is the main loop index, decreasing from N in steps of 2
*
         K = N
   10    CONTINUE
         KW = NB + K - N
*
*        Exit from loop
*
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LE.2 ) THEN
            IF ( NB.GE.N .AND. K.EQ.2 ) THEN
               CALL DCOPY( K-1, A( 1, K ), 1, W( 1, KW ), 1 )
               W( K, KW ) = ZERO
               IF( K.LT.N ) THEN
                  CALL DGEMV( 'No transpose', K, N-K, ONE,
     $                  A( 1, K+1 ), LDA, W( K, KW+1 ), LDW,
     $                  ONE, W( 1, KW ), 1 )
               END IF
               A( K-1, K ) = W( K-1, KW )
               IF ( ABS( A( K-1, K ) ) .EQ. ZERO) THEN
                  IF( INFO.EQ.0 )
     $               INFO = K
               END IF
               IPIV( K ) = 0
               K = K-2
            ELSEIF ( NB.GE.N .AND. K.EQ.1 ) THEN
               IF( INFO.EQ.0 )
     $            INFO = K
*               K = K-1
               KADJ = 1
            END IF
            GO TO 30
         END IF
*
*        Copy column K and K-1 of A to column K and K-1 of W and update them
*
         CALL DCOPY( K-1, A( 1, K ), 1, W( 1, KW ), 1 )
         CALL DCOPY( K-2, A( 1, K-1 ), 1, W( 1, KW-1 ), 1 )
         W( K, KW ) = ZERO
         W( K-1, KW-1 ) = ZERO
         IF( K.LT.N ) THEN
            CALL DGEMV( 'No transpose', K, N-K, ONE, A( 1, K+1 ),
     $              LDA, W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
            CALL DGEMV( 'No transpose', K-1, N-K, ONE, A( 1, K+1 ),
     $              LDA, W( K-1, KW+1 ), LDW, ONE, W( 1, KW-1 ), 1 )
         END IF

         W( K, KW-1 ) = -W( K-1, KW )
*
*        Determine rows and columns to be interchanged
*
         ABSAKP1K = ABS( W( K-1, KW ) )
*
*        IMAX1 is the row-index of the absolute value largest element in
*        row 1 to K-2, column K.
*        IMAX2 is the row-index of the absolute value largest element in
*        row 1 to K-2 column K-1.
*        COLMAX1 and COLMAX2 are their absolute values.
*
         IF(K.GT.2) THEN
            IMAX1 = IDAMAX( K-2, W( 1, KW ), 1 )
            COLMAX1 = ABS( W( IMAX1, KW ) )
            IMAX2 = IDAMAX( K-2, W( 1, KW-1 ), 1 )
            COLMAX2 = ABS( W( IMAX2, KW-1 ) )
         ELSE
            IMAX1 = 0
            COLMAX1 = ZERO
            IMAX2 = 0
            COLMAX2 = ZERO
         ENDIF
*
         IF( MAX(MAX( ABSAKP1K, COLMAX1 ), COLMAX2).EQ.ZERO ) THEN
*
*           Column K and K+1 is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = 0
            IPIV( K ) = KP
         ELSE
            IF( ABSAKP1K.GE.MAX( COLMAX1, COLMAX2 ) ) THEN
*
*              No interchange
*
               KP = 0
               IPIV( K ) = KP
            ELSE

               IF( COLMAX1.GE.COLMAX2 ) THEN

*
*                 Absolute value largest element is in column K
*                 Interchange rows and columns K-1 and IMAX1
*                  
                  KP = IMAX1
                  IPIV( K ) = KP

*                 
*                 Write the column KW-1 of W with elements in column IMAX1
*        
                  CALL DCOPY( IMAX1-1, A( 1, IMAX1 ), 1,
     $                     W( 1, KW-1 ), 1 )

                  W( IMAX1, KW-1 ) = ZERO

                  CALL DCOPY( K-IMAX1, A( IMAX1, IMAX1+1 ), LDA,
     $                     W( IMAX1+1, KW-1 ), 1 )

                  CALL DSCAL( K-IMAX1, -ONE, W( IMAX1+1, KW-1 ), 1)

*                 
*                 Update the column KW-1 of W
*     
                  IF( K.LT.N ) THEN
                     CALL DGEMV( 'No transpose', K, N-K, ONE,
     $                     A( 1, K+1 ), LDA, W( IMAX1, KW+1 ), LDW,
     $                     ONE, W( 1, KW-1 ), 1 )
                  END IF

*                  W( K, KW-1 ) = -W( K-1, KW )

*                 
*                 Write the column IMAX1 of A with elements in column K-1 of A
*
                  CALL DCOPY( IMAX1-1, A( 1, K-1 ), 1,
     $                     A( 1, IMAX1 ), 1 ) 

                  CALL DCOPY( K-IMAX1-2, A( IMAX1+1, K-1 ), 1,
     $                     A( IMAX1, IMAX1+1 ), LDA )

                  CALL DSCAL( K-IMAX1-2, -ONE, A( IMAX1, IMAX1+1 ),
     $                     LDA)
*                 
*                 Interchange rows K-1 and IMAX1 in last K-1 columns of A
*
                  IF( K.LT.N ) THEN
                     CALL DSWAP( N-K, A( K-1, K+1 ), LDA,
     $                        A( IMAX1, K+1 ), LDA )
                  END IF

*                 
*                 Interchange rows K-1 and IMAX1 in last KW-1 columns of W
*
                  CALL DSWAP( N-K+2, W( K-1, KW-1 ), LDW, 
     $                    W( IMAX1, KW-1 ), LDW )

               ELSE

*
*                 Absolute value largest element is in column K-1
*                 Interchange rows and columns K and K-1, then Interchange K-1 and IMAX2
*                                  
                  KP = -IMAX2
                  IPIV( K ) = KP

*                 
*                 Interchange columns KW and KW-1, then write the column KW-1 of W with elements in column IMAX2
*       
                  CALL DSWAP( K, W( 1, KW ), 1, W( 1, KW-1 ),
     $                     1 )

                  CALL DCOPY( IMAX2-1, A( 1, IMAX2 ), 1,
     $                     W( 1, KW-1 ), 1 )

                  W( IMAX2, KW-1 ) = ZERO

                  CALL DCOPY( K-IMAX2, A( IMAX2, IMAX2+1 ), LDA,
     $                     W( IMAX2+1, KW-1 ), 1 )

                  CALL DSCAL( K-IMAX2, -ONE, W( IMAX2+1, KW-1 ), 1)

*                 
*                 Update the column KW-1 of W
*     
                  IF( K.LT.N ) THEN
                     CALL DGEMV( 'No transpose', K, N-K, ONE,
     $                     A( 1, K+1 ), LDA, W( IMAX2, KW+1 ), LDW,
     $                     ONE, W( 1, KW-1 ), 1 )
                  END IF

*                  W( K, KW-1 ) = -W( K-1, KW )

*                 Interchange rows K and K-1 columns of A
*
                  CALL DSWAP( K-2, A( 1, K ), 1, A( 1, K-1 ),
     $                     1 )

                  A( K-1, K ) = -A( K-1, K )

*                 
*                 Write the column IMAX2 of A with elements in column K-1 of A
*
                  CALL DCOPY( IMAX2-1, A( 1, K-1 ), 1,
     $                     A( 1, IMAX2 ), 1 ) 

                  CALL DCOPY( K-IMAX2-2, A( IMAX2+1, K-1 ), 1,
     $                     A( IMAX2, IMAX2+1 ), LDA )

                  CALL DSCAL( K-IMAX2-2, -ONE, A( IMAX2, IMAX2+1 ),
     $					   LDA)
*                 
*                 Interchange rows K and K-1, then K-1 and IMAX2 in last K+1 columns of A
*
                  IF( K.LT.N ) THEN
                     CALL DSWAP( N-K, A( K, K+1 ), LDA,
     $                        A( K-1, K+1 ), LDA )

                     CALL DSWAP( N-K, A( K-1, K+1 ), LDA,
     $                        A( IMAX2, K+1 ), LDA )
                  END IF

*                 
*                 Interchange rows K and K-1, then K-1 and IMAX2 in last K-1 columns of W
*
                  CALL DSWAP( N-K+2, W( K, KW-1 ), LDW,
     $                     W( K-1, KW-1 ), LDW )

                  CALL DSWAP( N-K+2, W( K-1, KW-1 ), LDW,
     $                     W( IMAX2, KW-1 ), LDW )

               END IF
            END IF   

*                 
*           Write back C*S^-1 to A
*            
            DO 20 J = 1, K-2
               A( J, K-1 ) = W( J, KW )/W( K-1, KW )
               A( J, K ) = -W( J, KW-1 )/W( K-1, KW )
20             CONTINUE

            A( K-1, K ) = W( K-1, KW )

         END IF

         K = K-2

         GO TO 10
*
30    CONTINUE

         KW = NB + K - N
*
*        Update the upper triangle of A11 (= A(1:k,1:k)) as
*
*        A11 := A11 + U12*D*U12**T = A11 + U12*W**T
*
*        computing blocks of NB columns at a time
*
         DO 50 J = 1, K, NB
            JB = MIN( NB, K-J+1 )

*
*           Update the rectangular subdiagonal block
*
            IF( J+JB.LE.K )
     $         CALL DGEMM( 'No transpose', 'Transpose', K-J-JB+1,
     $                     JB, N-K, ONE, A( 1, K+1 ), LDA,
     $                     W( K-J-JB+2, KW+1 ), LDW, ONE,
     $                     A( 1, K-J-JB+2 ), LDA )
*
*           Update the upper triangle of the diagonal block
*
            DO 40 JJ = 1, JB - 1
               CALL DGEMV( 'No transpose', JJ, N-K, ONE,
     $                     A( K-J-JB+2, K+1 ), LDA,
     $                     W( K+JJ-J-JB+2, KW+1 ), LDW, ONE,
     $                     A( K-J-JB+2, K+JJ-J-JB+2 ), 1 )
  40       CONTINUE

  50    CONTINUE
*
*        Put U12 in standard form by partially undoing the interchanges
*        of rows in columns 1:k-1 looping backwards from k-1 to 1
*
         J = N - K - 1
  60    CONTINUE
*
*           Undo the interchanges (if any) of rows JJ and JP at each
*           step J
*
*           (Here, J is a diagonal index)

            IF( J.GT.1 ) THEN
               JJ = N-J+1
               JP = IPIV( N-J+1 )

               IF( JP.LT.0 ) THEN
                  JP = -JP
*                 (Here, J is a diagonal index)
                  CALL DSWAP( J-1, A( JP, N-J+2 ), LDA,
     $                  A( JJ-1, N-J+2 ), LDA )
                  CALL DSWAP( J-1, A( JJ-1, N-J+2 ), LDA,
     $                  A( JJ, N-J+2 ), LDA )
               ELSEIF( JP.GT.0 ) THEN
                  CALL DSWAP( J-1, A( JP, N-J+2 ), LDA,
     $                  A( JJ-1, N-J+2 ), LDA )
               END IF
               
            END IF
*        (NOTE: Here, J is used to determine row length. Length J
*        of the rows to swap back doesn't include diagonal element)

         J = J - 2
         IF( J.GT.1 )
     $      GO TO 60
*
*        Set KB to the number of columns factorized
*
         KB = N - K + KADJ
*         
      ELSE
*
*        Factorize the leading columns of A using the lower triangle
*        of A and working forwards, and compute the matrix W = L21*D
*        for use in updating A22
*
*        K is the main loop index, increasing from 1 in steps 2
*
         K = 1
   70    CONTINUE
*
*        Exit from loop
*
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GE.N-1 ) THEN
            IF( NB.GE.N .AND. K.EQ.N-1 ) THEN
               CALL DCOPY( N-K, A( K+1, K ), 1, W( K+1, K ), 1 )
               W( K, K ) = ZERO
               CALL DGEMV( 'No transpose', N-K+1, K-1, ONE,
     $               A( K, 1 ), LDA, W( K, 1 ), LDW, ONE,
     $               W( K, K ), 1 )
               A( K+1, K ) = W( K+1, K )
               IF ( ABS( A( K+1, K ) ) .EQ. ZERO) THEN
                  IF( INFO.EQ.0 )
     $               INFO = K
               END IF
               IPIV( K ) = 0
               K = K+2
            ELSEIF( NB.GE.N .AND. K.EQ.N ) THEN
               IF( INFO.EQ.0 )
     $            INFO = K
*               K = K+1
               KADJ = 1
            END IF
            GO TO 90
         END IF
*
*        Copy column K and K+1 of A to column K and K+1 of W and update them
*
         CALL DCOPY( N-K, A( K+1, K ), 1, W( K+1, K ), 1 )
         CALL DCOPY( N-K-1, A( K+2, K+1 ), 1, W( K+2, K+1 ), 1 )
         W( K, K ) = ZERO
         W( K+1, K+1 ) = ZERO
         CALL DGEMV( 'No transpose', N-K+1, K-1, ONE, A( K, 1 ),
     $               LDA, W( K, 1 ), LDW, ONE, W( K, K ), 1 )
         CALL DGEMV( 'No transpose', N-K, K-1, ONE, A( K+1, 1 ),
     $               LDA, W( K+1, 1 ), LDW, ONE, W( K+1, K+1 ), 1 )

         W( K, K+1 ) = -W( K+1, K )
*
*        Determine rows and columns to be interchanged
*
         ABSAKP1K = ABS( W( K+1, K ) )
*
*        IMAX1 is the row-index of the absolute value largest element in
*        row K+2 to N, column K.
*        IMAX2 is the row-index of the absolute value largest element in
*        row K+2 to N, column K+1.
*        COLMAX1 and COLMAX2 are their absolute values.
*
         IF(K.LT.N-1) THEN
            IMAX1 = K+1 + IDAMAX( N-K-1, W( K+2, K ), 1 )
            COLMAX1 = ABS( W( IMAX1, K ) )
            IMAX2 = K+1 + IDAMAX( N-K-1, W( K+2, K+1 ), 1 )
            COLMAX2 = ABS( W( IMAX2, K+1 ) )
         ELSE
            IMAX1 = 0
            COLMAX1 = ZERO
            IMAX2 = 0
            COLMAX2 = ZERO
         ENDIF
*
         IF( MAX(MAX( ABSAKP1K, COLMAX1 ), COLMAX2).EQ.ZERO ) THEN
*
*           Column K and K+1 is zero or underflow: set INFO and continue
*
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = 0
            IPIV( K ) = KP
         ELSE
            IF( ABSAKP1K.GE.MAX( COLMAX1, COLMAX2 ) ) THEN
*
*              No interchange
*
               KP = 0
               IPIV( K ) = KP
            ELSE

               IF( COLMAX1.GE.COLMAX2 ) THEN

*
*                 Absolute value largest element is in column K
*                 Interchange rows and columns K+1 and IMAX1
*                  
                  KP = IMAX1
                  IPIV( K ) = KP

*                 
*                 Write the column K+1 of W with elements in column IMAX1
*        
                  CALL DCOPY( IMAX1-K, A( IMAX1, K ), LDA,
     $                     W( K, K+1 ), 1 )

                  CALL DSCAL( IMAX1-K, -ONE, W( K, K+1 ), 1)

                  W( IMAX1, K+1 ) = ZERO

                  CALL DCOPY( N-IMAX1, A( IMAX1+1, IMAX1 ), 1,
     $                     W( IMAX1+1, K+1 ), 1 )

*                 
*                 Update the column K+1 of W
*     
                  CALL DGEMV( 'No transpose', N-K+1, K-1, ONE,
     $                     A( K, 1 ), LDA, W( IMAX1, 1 ), LDW, ONE,
     $                     W( K, K+1 ), 1 )

*                  W( K, K+1 ) = -W( K+1, K )

*                 
*                 Write the column IMAX1 of A with elements in column K+1 of A
*
                  CALL DCOPY( IMAX1-K-2, A( K+2, K+1 ), 1,
     $                     A( IMAX1, K+2 ), LDA )

                  CALL DSCAL( IMAX1-K-2, -ONE, A( IMAX1, K+2 ), LDA)

                  CALL DCOPY( N-IMAX1, A( IMAX1+1, K+1 ), 1,
     $                     A( IMAX1+1, IMAX1 ), 1 ) 

*                 
*                 Interchange rows K+1 and IMAX1 in first K-1 columns of A
*
                  CALL DSWAP( K-1, A( K+1, 1 ), LDA, A( IMAX1, 1 ),
     $                     LDA )

*                 
*                 Interchange rows K+1 and IMAX1 in first K-1 columns of W
*
                  CALL DSWAP( K+1, W( K+1, 1 ), LDW, W( IMAX1, 1 ),
     $                     LDW )

               ELSE

*
*                 Absolute value largest element is in column K+1
*                 Interchange rows and columns K and K+1, then Interchange K+1 and IMAX2
*                                  
                  KP = -IMAX2
                  IPIV( K ) = KP

*                 
*                 Interchange columns K and K+1, then write the column K+1 of W with elements in column IMAX2
*       
                  CALL DSWAP( N-K+1, W( K, K ), 1, W( K, K+1 ),
     $                     1 )

                  CALL DCOPY( IMAX2-K, A( IMAX2, K ), LDA,
     $                     W( K, K+1 ), 1 )

                  CALL DSCAL( IMAX2-K, -ONE, W( K, K+1 ), 1)

                  W( IMAX2, K+1 ) = ZERO

                  CALL DCOPY( N-IMAX2, A( IMAX2+1, IMAX2 ), 1,
     $                     W( IMAX2+1, K+1 ), 1 ) 

*                 
*                 Update the column K+1 of W
*     
                  CALL DGEMV( 'No transpose', N-K+1, K-1, ONE,
     $                     A( K, 1 ), LDA, W( IMAX2, 1 ), LDW, ONE,
     $                     W( K, K+1 ), 1 )

*                  W( K, K+1 ) = -W( K+1, K )

*                 Interchange rows K and K+1 columns of A
*
                  CALL DSWAP( N-K-1, A( K+2, K ), 1, A( K+2, K+1 ),
     $                     1 )

                  A( K+1, K ) = -A( K+1, K )

*                 
*                 Write the column IMAX2 of A with elements in column K+1 of A
*
                  CALL DCOPY( IMAX2-K-2, A( K+2, K+1 ), 1,
     $                     A( IMAX2, K+2 ), LDA )

                  CALL DSCAL( IMAX2-K-2, -ONE, A( IMAX2, K+2 ), LDA)

                  CALL DCOPY( N-IMAX2, A( IMAX2+1, K+1 ), 1,
     $                     A( IMAX2+1, IMAX2 ), 1 ) 

*                 
*                 Interchange rows K and K+1, then K+1 and IMAX2 in first K-1 columns of A
*
                  CALL DSWAP( K-1, A( K, 1 ), LDA, A( K+1, 1 ),
     $                     LDA )

                  CALL DSWAP( K-1, A( K+1, 1 ), LDA, A( IMAX2, 1 ),
     $                     LDA )

*                 
*                 Interchange rows K and K+1, then K+1 and IMAX2 in first K-1 columns of W
*
                  CALL DSWAP( K+1, W( K, 1 ), LDW, W( K+1, 1 ),
     $                     LDW )

                  CALL DSWAP( K+1, W( K+1, 1 ), LDW, W( IMAX2, 1 ),
     $                     LDW )

               END IF
            END IF   

*                 
*           Write back C*S^-1 to A
*            
            DO 80 J = K+2, N
               A( J, K ) = -W( J, K+1 )/W( K+1, K )
               A( J, K+1 ) = W( J, K )/W( K+1, K )
80             CONTINUE

            A( K+1, K ) = W( K+1, K )

         END IF

         K = K+2

         GO TO 70
*
90    CONTINUE
*
*        Update the lower triangle of A22 (= A(k:n,k:n)) as
*
*        A22 := A22 + L21*D*L21**T = A22 + L21*W**T
*
*        computing blocks of NB columns at a time
*
         DO 110 J = K, N, NB
            JB = MIN( NB, N-J+1 )
*
*           Update the lower triangle of the diagonal block
*
            DO 100 JJ = J, J + JB - 2
               CALL DGEMV( 'No transpose', J+JB-JJ-1, K-1, ONE,
     $                     A( JJ+1, 1 ), LDA, W( JJ, 1 ), LDW,
     $                     ONE, A( JJ+1, JJ ), 1 )
  100       CONTINUE
*
*           Update the rectangular subdiagonal block
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'Transpose', N-J-JB+1,
     $                     JB, K-1, ONE, A( J+JB, 1 ), LDA,
     $                     W( J, 1 ), LDW, ONE, A( J+JB, J ),
     $                     LDA )
  110    CONTINUE
*
*        Put L21 in standard form by partially undoing the interchanges
*        of rows in columns 1:k-1 looping backwards from k-1 to 1
*
         J = K - 2
  120    CONTINUE
*
*           Undo the interchanges (if any) of rows JJ and JP at each
*           step J
*
*           (Here, J is a diagonal index)

            IF( J.GT.1 ) THEN
               JJ = J
               JP = IPIV( J )

               IF( JP.LT.0 ) THEN
                  JP = -JP
*                 (Here, J is a diagonal index)
                  CALL DSWAP( J-1, A( JP, 1 ), LDA, A( JJ+1, 1 ),
     $				   LDA )
                  CALL DSWAP( J-1, A( JJ+1, 1 ), LDA, A( JJ, 1 ),
     $				   LDA )
               ELSEIF( JP.GT.0 ) THEN
                  CALL DSWAP( J-1, A( JP, 1 ), LDA, A( JJ+1, 1 ),
     $				   LDA )
               END IF
               
            END IF
*        (NOTE: Here, J is used to determine row length. Length J
*        of the rows to swap back doesn't include diagonal element)

         J = J - 2
         IF( J.GT.1 )
     $      GO TO 120
*
*        Set KB to the number of columns factorized
*
         KB = K - 1 + KADJ
*
      END IF
      RETURN
*
*     End of SLASYF
*
      END
