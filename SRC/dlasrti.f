*> \brief \b DLASRTI sorts array indices based on the referenced numbers
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLASRTI + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlasrti.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlasrti.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlasrti.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLASRTI( ID, N, X, INDICES, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          ID
*       INTEGER            INFO, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   X( * )
*       INTEGER            INDICES( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> Sort the numbers in X indirectly in increasing order (if ID = 'I') or
*> in decreasing order (if ID = 'D' ) using the array of INDICES.
*>
*> Use Quick Sort, reverting to Insertion sort on arrays of
*> size <= 20. Dimension of STACK limits N to about 2**32.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] ID
*> \verbatim
*>          ID is CHARACTER*1
*>          = 'I': sort X in increasing order;
*>          = 'D': sort X in decreasing order.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The length of the array X and INDICES.
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array, dimension (N)
*>          The values to be sorted.
*> \endverbatim
*>
*> \param[in,out] INDICES
*> \verbatim
*>          INDICES is INTEGER array, dimension (N)
*>          On entry, the indices of values in X to be sorted.
*>          On exit, the indices have been sorted such that
*>          X( INDICES(1) ) <= ... <= X( INDICES(N) )
*>          or decreasing order such that
*>          X( INDICES(1) ) >= ... >= X( INDICES(N) )
*>          depending on ID.
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
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup auxOTHERcomputational
*
*  =====================================================================
      SUBROUTINE DLASRTI( ID, N, X, INDICES, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          ID
      INTEGER            INFO, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   X( * )
      INTEGER            INDICES( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            SELECT
      PARAMETER          ( SELECT = 20 )
*     ..
*     .. Local Scalars ..
      INTEGER            DIR, ENDD, I, J, START, STKPNT, TMP
      DOUBLE PRECISION   P1, P2, P3, PIVOT
*     ..
*     .. Local Arrays ..
      INTEGER            STACK( 2, 32 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      DIR = -1
      IF( LSAME( ID, 'D' ) ) THEN
         DIR = 0
      ELSE IF( LSAME( ID, 'I' ) ) THEN
         DIR = 1
      END IF
      IF( DIR.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DLASRTI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
      STKPNT = 1
      STACK( 1, 1 ) = 1
      STACK( 2, 1 ) = N
   10 CONTINUE
      START = STACK( 1, STKPNT )
      ENDD = STACK( 2, STKPNT )
      STKPNT = STKPNT - 1
      IF( ENDD-START.LE.SELECT .AND. ENDD-START.GT.0 ) THEN
*
*        Do Insertion sort on X( START:ENDD )
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            DO 30 I = START + 1, ENDD
               DO 20 J = I, START + 1, -1
                  IF( X( INDICES(J) ).GT.X( INDICES(J-1) ) ) THEN
                     TMP = INDICES( J )
                     INDICES( J ) = INDICES( J-1 )
                     INDICES( J-1 ) = TMP
                  ELSE
                     GO TO 30
                  END IF
   20          CONTINUE
   30       CONTINUE
*
         ELSE
*
*           Sort into increasing order
*
            DO 50 I = START + 1, ENDD
               DO 40 J = I, START + 1, -1
                  IF( X( INDICES(J) ).LT.X( INDICES(J-1) ) ) THEN
                     TMP = INDICES( J )
                     INDICES( J ) = INDICES( J-1 )
                     INDICES( J-1 ) = TMP
                  ELSE
                     GO TO 50
                  END IF
   40          CONTINUE
   50       CONTINUE
*
         END IF
*
      ELSE IF( ENDD-START.GT.SELECT ) THEN
*
*        Partition X( START:ENDD ) and stack parts, largest one first
*
*        Choose partition entry as median of 3
*
         P1 = X( INDICES(START) )
         P2 = X( INDICES(ENDD) )
         I = ( START+ENDD ) / 2
         P3 = X( INDICES(I) )
         IF( P1.LT.P2 ) THEN
            IF( P3.LT.P1 ) THEN
               PIVOT = P1
            ELSE IF( P3.LT.P2 ) THEN
               PIVOT = P3
            ELSE
               PIVOT = P2
            END IF
         ELSE
            IF( P3.LT.P2 ) THEN
               PIVOT = P2
            ELSE IF( P3.LT.P1 ) THEN
               PIVOT = P3
            ELSE
               PIVOT = P1
            END IF
         END IF
*
         IF( DIR.EQ.0 ) THEN
*
*           Sort into decreasing order
*
            I = START - 1
            J = ENDD + 1
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF( X( INDICES(J) ).LT.PIVOT )
     $         GO TO 70
   80       CONTINUE
            I = I + 1
            IF( X( INDICES(I) ).GT.PIVOT )
     $         GO TO 80
            IF( I.LT.J ) THEN
               TMP = INDICES( I )
               INDICES( I ) = INDICES( J )
               INDICES( J ) = TMP
               GO TO 60
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         ELSE
*
*           Sort into increasing order
*
            I = START - 1
            J = ENDD + 1
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF( X( INDICES(J) ).GT.PIVOT )
     $         GO TO 100
  110       CONTINUE
            I = I + 1
            IF( X( INDICES(I) ).LT.PIVOT )
     $         GO TO 110
            IF( I.LT.J ) THEN
               TMP = INDICES( I )
               INDICES( I ) = INDICES( J )
               INDICES( J ) = TMP
               GO TO 90
            END IF
            IF( J-START.GT.ENDD-J-1 ) THEN
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
            ELSE
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = J + 1
               STACK( 2, STKPNT ) = ENDD
               STKPNT = STKPNT + 1
               STACK( 1, STKPNT ) = START
               STACK( 2, STKPNT ) = J
            END IF
         END IF
      END IF
      IF( STKPNT.GT.0 )
     $   GO TO 10
      RETURN
*
*     End of DLASRTI
*
      END
