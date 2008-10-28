      INTEGER FUNCTION ILA_LEN_TRIM(SUBNAM)
C
C  -- LAPACK auxiliary routine (version 3.1) --
C     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
C     October 2006
C
C     .. Scalar Arguments ..
      CHARACTER*(*) SUBNAM
C     ..
C
C  Purpose
C  =======
C
C  ILA_LEN_TRIM is called from testing and timing routines to remove
C  trailing spaces from its argument.  It is included in the library
C  for possible use within a user's XERBLA error-handing routine.
C
C  Arguments
C  =========
C
C  SUBNAM  (input) CHARACTER*(*)
C          Provides the string.
C
C  RETURN VALUE:  INTEGER
C          = N > 0 : The location of the last non-blank.
C          = 0     : The entire string is blank.
C
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC LEN
C     ..

      DO I = LEN(SUBNAM),1,-1
          IF (SUBNAM(I:I).NE.' ') THEN
              ILA_LEN_TRIM = I
              RETURN
          END IF
      END DO
      ILA_LEN_TRIM = 0
      END
