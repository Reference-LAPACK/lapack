*> \brief \b SLAHRD
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition
*  ==========
*
*       SUBROUTINE SLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
* 
*       .. Scalar Arguments ..
*       INTEGER            K, LDA, LDT, LDY, N, NB
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), T( LDT, NB ), TAU( NB ),
*      $                   Y( LDY, NB )
*       ..
*  
*  Purpose
*  =======
*
*>\details \b Purpose:
*>\verbatim
*>
*> SLAHRD reduces the first NB columns of a real general n-by-(n-k+1)
*> matrix A so that elements below the k-th subdiagonal are zero. The
*> reduction is performed by an orthogonal similarity transformation
*> Q**T * A * Q. The routine returns the matrices V and T which determine
*> Q as a block reflector I - V*T*V**T, and also the matrix Y = A * V * T.
*>
*> This is an OBSOLETE auxiliary routine. 
*> This routine will be 'deprecated' in a  future release.
*> Please use the new routine SLAHR2 instead.
*>
*>\endverbatim
*
*  Arguments
*  =========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The offset for the reduction. Elements below the k-th
*>          subdiagonal in the first NB columns are reduced to zero.
*> \endverbatim
*>
*> \param[in] NB
*> \verbatim
*>          NB is INTEGER
*>          The number of columns to be reduced.
*> \endverbatim
*>
*
*  Authors
*  =======
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup realOTHERauxiliary
*
*
*  Further Details
*  ===============
*>\details \b Further \b Details
*> \verbatim
*          unchanged. See Further Details.
*>
*>  LDA     (input) INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*>
*>  TAU     (output) REAL array, dimension (NB)
*>          The scalar factors of the elementary reflectors. See Further
*>          Details.
*>
*>  T       (output) REAL array, dimension (LDT,NB)
*>          The upper triangular matrix T.
*>
*>  LDT     (input) INTEGER
*>          The leading dimension of the array T.  LDT >= NB.
*>
*>  Y       (output) REAL array, dimension (LDY,NB)
*>          The n-by-nb matrix Y.
*>
*>  LDY     (input) INTEGER
*>          The leading dimension of the array Y. LDY >= N.
*>
*>
*>  The matrix Q is represented as a product of nb elementary reflectors
*>
*>     Q = H(1) H(2) . . . H(nb).
*>
*>  Each H(i) has the form
*>
*>     H(i) = I - tau * v * v**T
*>
*>  where tau is a real scalar, and v is a real vector with
*>  v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
*>  A(i+k+1:n,i), and tau in TAU(i).
*>
*>  The elements of the vectors v together form the (n-k+1)-by-nb matrix
*>  V which is needed, with T and Y, to apply the transformation to the
*>  unreduced part of the matrix, using an update of the form:
*>  A := (I - V*T*V**T) * (A - Y*V**T).
*>
*>  The contents of A on exit are illustrated by the following example
*>  with n = 7, k = 3 and nb = 2:
*>
*>     ( a   h   a   a   a )
*>     ( a   h   a   a   a )
*>     ( a   h   a   a   a )
*>     ( h   h   a   a   a )
*>     ( v1  h   a   a   a )
*>     ( v1  v2  a   a   a )
*>     ( v1  v2  a   a   a )
*>
*>  where a denotes an element of the original matrix A, h denotes a
*>  modified element of the upper Hessenberg matrix H, and vi denotes an
*>  element of the vector defining H(i).
*>
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE SLAHRD( N, K, NB, A, LDA, TAU, T, LDT, Y, LDY )
*
*  -- LAPACK auxiliary routine (version 3.3.1) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER            K, LDA, LDT, LDY, N, NB
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), T( LDT, NB ), TAU( NB ),
     $                   Y( LDY, NB )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               EI
*     ..
*     .. External Subroutines ..
      EXTERNAL           SAXPY, SCOPY, SGEMV, SLARFG, SSCAL, STRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      DO 10 I = 1, NB
         IF( I.GT.1 ) THEN
*
*           Update A(1:n,i)
*
*           Compute i-th column of A - Y * V**T
*
            CALL SGEMV( 'No transpose', N, I-1, -ONE, Y, LDY,
     $                  A( K+I-1, 1 ), LDA, ONE, A( 1, I ), 1 )
*
*           Apply I - V * T**T * V**T to this column (call it b) from the
*           left, using the last column of T as workspace
*
*           Let  V = ( V1 )   and   b = ( b1 )   (first I-1 rows)
*                    ( V2 )             ( b2 )
*
*           where V1 is unit lower triangular
*
*           w := V1**T * b1
*
            CALL SCOPY( I-1, A( K+1, I ), 1, T( 1, NB ), 1 )
            CALL STRMV( 'Lower', 'Transpose', 'Unit', I-1, A( K+1, 1 ),
     $                  LDA, T( 1, NB ), 1 )
*
*           w := w + V2**T *b2
*
            CALL SGEMV( 'Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ),
     $                  LDA, A( K+I, I ), 1, ONE, T( 1, NB ), 1 )
*
*           w := T**T *w
*
            CALL STRMV( 'Upper', 'Transpose', 'Non-unit', I-1, T, LDT,
     $                  T( 1, NB ), 1 )
*
*           b2 := b2 - V2*w
*
            CALL SGEMV( 'No transpose', N-K-I+1, I-1, -ONE, A( K+I, 1 ),
     $                  LDA, T( 1, NB ), 1, ONE, A( K+I, I ), 1 )
*
*           b1 := b1 - V1*w
*
            CALL STRMV( 'Lower', 'No transpose', 'Unit', I-1,
     $                  A( K+1, 1 ), LDA, T( 1, NB ), 1 )
            CALL SAXPY( I-1, -ONE, T( 1, NB ), 1, A( K+1, I ), 1 )
*
            A( K+I-1, I-1 ) = EI
         END IF
*
*        Generate the elementary reflector H(i) to annihilate
*        A(k+i+1:n,i)
*
         CALL SLARFG( N-K-I+1, A( K+I, I ), A( MIN( K+I+1, N ), I ), 1,
     $                TAU( I ) )
         EI = A( K+I, I )
         A( K+I, I ) = ONE
*
*        Compute  Y(1:n,i)
*
         CALL SGEMV( 'No transpose', N, N-K-I+1, ONE, A( 1, I+1 ), LDA,
     $               A( K+I, I ), 1, ZERO, Y( 1, I ), 1 )
         CALL SGEMV( 'Transpose', N-K-I+1, I-1, ONE, A( K+I, 1 ), LDA,
     $               A( K+I, I ), 1, ZERO, T( 1, I ), 1 )
         CALL SGEMV( 'No transpose', N, I-1, -ONE, Y, LDY, T( 1, I ), 1,
     $               ONE, Y( 1, I ), 1 )
         CALL SSCAL( N, TAU( I ), Y( 1, I ), 1 )
*
*        Compute T(1:i,i)
*
         CALL SSCAL( I-1, -TAU( I ), T( 1, I ), 1 )
         CALL STRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T, LDT,
     $               T( 1, I ), 1 )
         T( I, I ) = TAU( I )
*
   10 CONTINUE
      A( K+NB, NB ) = EI
*
      RETURN
*
*     End of SLAHRD
*
      END
