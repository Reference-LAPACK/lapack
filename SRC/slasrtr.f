*> \brief \b SLASRTR sorts array indices based on the referenced numbers
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SLASRTR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasrtr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasrtr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasrtr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SLASRTR( ID, M, N, A, LDA, IPVT, WORK, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          ID
*       INTEGER            INFO, M, N, LDA
*       ..
*       .. Array Arguments ..
*       REAL               A( LDA, * ), WORK( * )
*       INTEGER            IPVT( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*> Apply row sorting based on the maximum norm of the rows of A. The
*> rows can be sorted in increasing (if ID = 'I') or decreasing order
*> (if ID = 'D').
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ID
*> \verbatim
*>          ID is CHARACTER*1
*>          = 'I': sort the rows of A by increasing row norm;
*>          = 'D': sort the rows of A by decreasing row norm.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A. N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          On entry, the M by N matrix A.
*>          On exit, A contains the row-permuted matrix A.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[out] IPVT
*> \verbatim
*>          IPVT is INTEGER array, dimension (M)
*>          On exit, if IPVT(J)=K, then the J-th row of A*P was the
*>          the K-th row of A.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is REAL array, dimension (M).
*>          On exit, WORK contains the maximum norms of the rows of A.
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
*> \author Christoph Conrads (https://christoph-conrads.name)
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
      SUBROUTINE SLASRTR( ID, M, N, A, LDA, IPVT, WORK, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, M, N, LDA
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * ), WORK( * )
      INTEGER            IPVT( * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, SLAPMR, SLASRTI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( .NOT.LSAME( ID, 'I' ) .AND. .NOT.LSAME( ID, 'D' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SLASRTR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.1 .OR. N.LT.1 )
     $   RETURN
*
*     Compute maximum norm of each row
*
      DO I = 1, M
         WORK( I ) = ABS( A( I, 1 ) )
      END DO
*
      DO I = 1, M
         DO J = 2, N
            WORK( I ) = MAX( WORK( I ), ABS( A( I, J ) ) )
         END DO
      END DO
*
*     Sort row indices
*
      DO I = 1, M
         IPVT( I ) = I
      END DO
*
      CALL SLASRTI( ID, M, WORK, IPVT, INFO )
      IF( INFO.NE.0 ) THEN
         RETURN
      END IF
*
*     Sort rows
*
      CALL SLAPMR( .TRUE., M, N, A, LDA, IPVT )
*
*     End of SLASRTR
*
      END
