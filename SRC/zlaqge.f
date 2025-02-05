*> \brief \b ZLAQGE scales a general rectangular matrix, using row and column scaling factors computed by sgeequ.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download ZLAQGE + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqge.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqge.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqge.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
*                          EQUED )
*
*       .. Scalar Arguments ..
*       CHARACTER          EQUED
*       INTEGER            LDA, M, N
*       DOUBLE PRECISION   AMAX, COLCND, ROWCND
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   C( * ), R( * )
*       COMPLEX*16         A( LDA, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLAQGE equilibrates a general M by N matrix A using the row and
*> column scaling factors in the vectors R and C.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is COMPLEX*16 array, dimension (LDA,N)
*>          On entry, the M by N matrix A.
*>          On exit, the equilibrated matrix.  See EQUED for the form of
*>          the equilibrated matrix.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(M,1).
*> \endverbatim
*>
*> \param[in] R
*> \verbatim
*>          R is DOUBLE PRECISION array, dimension (M)
*>          The row scale factors for A.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (N)
*>          The column scale factors for A.
*> \endverbatim
*>
*> \param[in] ROWCND
*> \verbatim
*>          ROWCND is DOUBLE PRECISION
*>          Ratio of the smallest R(i) to the largest R(i).
*> \endverbatim
*>
*> \param[in] COLCND
*> \verbatim
*>          COLCND is DOUBLE PRECISION
*>          Ratio of the smallest C(i) to the largest C(i).
*> \endverbatim
*>
*> \param[in] AMAX
*> \verbatim
*>          AMAX is DOUBLE PRECISION
*>          Absolute value of largest matrix entry.
*> \endverbatim
*>
*> \param[out] EQUED
*> \verbatim
*>          EQUED is CHARACTER*1
*>          Specifies the form of equilibration that was done.
*>          = 'N':  No equilibration
*>          = 'R':  Row equilibration, i.e., A has been premultiplied by
*>                  diag(R).
*>          = 'C':  Column equilibration, i.e., A has been postmultiplied
*>                  by diag(C).
*>          = 'B':  Both row and column equilibration, i.e., A has been
*>                  replaced by diag(R) * A * diag(C).
*> \endverbatim
*
*> \par Internal Parameters:
*  =========================
*>
*> \verbatim
*>  THRESH is a threshold value used to decide if row or column scaling
*>  should be done based on the ratio of the row or column scaling
*>  factors.  If ROWCND < THRESH, row scaling is done, and if
*>  COLCND < THRESH, column scaling is done.
*>
*>  LARGE and SMALL are threshold values used to decide if row scaling
*>  should be done based on the absolute size of the largest matrix
*>  element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
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
*> \ingroup laqge
*
*  =====================================================================
      SUBROUTINE ZLAQGE( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX,
     $                   EQUED )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          EQUED
      INTEGER            LDA, M, N
      DOUBLE PRECISION   AMAX, COLCND, ROWCND
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), R( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, THRESH
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   CJ, LARGE, SMALL
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF
*
*     Initialize LARGE and SMALL.
*
      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL
*
      IF( ROWCND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE )
     $     THEN
*
*        No row scaling
*
         IF( COLCND.GE.THRESH ) THEN
*
*           No column scaling
*
            EQUED = 'N'
         ELSE
*
*           Column scaling
*
            DO 20 J = 1, N
               CJ = C( J )
               DO 10 I = 1, M
                  A( I, J ) = CJ*A( I, J )
   10          CONTINUE
   20       CONTINUE
            EQUED = 'C'
         END IF
      ELSE IF( COLCND.GE.THRESH ) THEN
*
*        Row scaling, no column scaling
*
         DO 40 J = 1, N
            DO 30 I = 1, M
               A( I, J ) = R( I )*A( I, J )
   30       CONTINUE
   40    CONTINUE
         EQUED = 'R'
      ELSE
*
*        Row and column scaling
*
         DO 60 J = 1, N
            CJ = C( J )
            DO 50 I = 1, M
               A( I, J ) = CJ*R( I )*A( I, J )
   50       CONTINUE
   60    CONTINUE
         EQUED = 'B'
      END IF
*
      RETURN
*
*     End of ZLAQGE
*
      END
