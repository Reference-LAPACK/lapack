      SUBROUTINE CLAG2Z( M, N, SA, LDSA, A, LDA, INFO)
*
*  -- LAPACK PROTOTYPE auxilary routine (version 3.1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     January 2007
*
*     ..
*     .. WARNING: PROTOTYPE ..
*     This is an LAPACK PROTOTYPE routine which means that the
*     interface of this routine is likely to be changed in the future
*     based on community feedback.
*
*     ..
*     .. Scalar Arguments ..
      INTEGER INFO,LDA,LDSA,M,N
*     ..
*     .. Array Arguments ..
      COMPLEX SA(LDSA,*)
      COMPLEX*16 A(LDA,*)
*     ..
*
*  Purpose
*  =======
*
*  CLAG2Z converts a COMPLEX SINGLE PRECISION matrix, SA, to a COMPLEX
*  DOUBLE PRECISION matrix, A.
*
*  Note that while it is possible to overflow while converting 
*  from double to single, it is not possible to overflow when
*  converting from single to double. 
*
*  This is a helper routine so there is no argument checking.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of lines of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  SA      (output) REAL array, dimension (LDSA,N)
*          On exit, the M-by-N coefficient matrix SA.
*
*  LDSA    (input) INTEGER
*          The leading dimension of the array SA.  LDSA >= max(1,M).
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N coefficient matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*  =========
*
*     .. Local Scalars ..
      INTEGER I,J
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      DO 20 J = 1,N
          DO 30 I = 1,M
              A(I,J) = SA(I,J)
   30     CONTINUE
   20 CONTINUE
      RETURN
*
*     End of CLAG2Z
*
      END
