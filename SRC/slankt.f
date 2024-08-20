*> \brief \b SLANKT returns the value of the 1-norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a real skew-symmetric tridiagonal matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SLANKT + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slankt.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slankt.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slankt.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       REAL             FUNCTION SLANKT( NORM, N, E )
*
*       .. Scalar Arguments ..
*       CHARACTER          NORM
*       INTEGER            N
*       ..
*       .. Array Arguments ..
*       REAL               E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SLANKT  returns the value of the one norm,  or the Frobenius norm, or
*> the  infinity norm,  or the  element of  largest absolute value  of a
*> real skew-symmetric tridiagonal matrix A.
*> \endverbatim
*>
*> \return SLANKT
*> \verbatim
*>
*>    SLANKT = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*>             (
*>             ( norm1(A),         NORM = '1', 'O' or 'o'
*>             (
*>             ( normI(A),         NORM = 'I' or 'i'
*>             (
*>             ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*>
*> where  norm1  denotes the  one norm of a matrix (maximum column sum),
*> normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*> normF  denotes the  Frobenius norm of a matrix (square root of sum of
*> squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] NORM
*> \verbatim
*>          NORM is CHARACTER*1
*>          Specifies the value to be returned in SLANKT as described
*>          above.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.  When N = 0, SLANKT is
*>          set to zero.
*> \endverbatim
*>
*> \param[in] E
*> \verbatim
*>          E is REAL array, dimension (N-1)
*>          The (n-1) sub-diagonal or super-diagonal elements of A.
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
*> \ingroup lankt
*
*  =====================================================================
      REAL             FUNCTION SLANKT( NORM, N, E )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            N
*     ..
*     .. Array Arguments ..
      REAL               E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      REAL               ANORM, SCALE, SUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME, SISNAN
      EXTERNAL           LSAME, SISNAN
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         ANORM = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         ANORM = ABS( E( N-1 ) )
         DO 10 I = 1, N - 2
            SUM = ABS( E( I ) )
            IF( ANORM .LT. SUM .OR. SISNAN( SUM ) ) ANORM = SUM
   10    CONTINUE
      ELSE IF( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' .OR.
     $         LSAME( NORM, 'I' ) ) THEN
*
*        Find norm1(A).
*
         IF( N.EQ.1 ) THEN
            ANORM = ZERO
         ELSE
            ANORM = ABS( E( 1 ) )
            SUM = ABS( E( N-1 ) )
            IF( ANORM .LT. SUM .OR. SISNAN( SUM ) ) ANORM = SUM
            DO 20 I = 2, N - 1
               SUM = ABS( E( I ) )+ABS( E( I-1 ) )
               IF( ANORM .LT. SUM .OR. SISNAN( SUM ) ) ANORM = SUM
   20       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR.
     $         ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( N.GT.1 ) THEN
            CALL SLASSQ( N-1, E, 1, SCALE, SUM )
            SUM = 2*SUM
         END IF
         ANORM = SCALE*SQRT( SUM )
      END IF
*
      SLANKT = ANORM
      RETURN
*
*     End of SLANKT
*
      END
