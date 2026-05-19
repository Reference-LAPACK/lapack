*> \brief \b SLAREV
*
* =========== DOCUMENTATION ===========
*
*> \par Purpose:
* =============
*>
*> \details \b Purpose:
*> \verbatim
*>
*> SLAREV reverses the order of the rows and/or columns of a given
*> M-by-N single precision matrix A.
*> \endverbatim
*
* Arguments:
* ==========
*
*> \param[in] WHICH
*> \verbatim
*>          WHICH is CHARACTER*1
*>          Specifies whether to reverse rows, columns, or both:
*>          = 'B': Reverse both rows and columns.
*>          = 'R': Reverse rows only.
*>          = 'C': Reverse columns only.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix A.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix A.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA, N)
*>          On entry, the matrix A to be reversed.
*>          On exit, the elements of A have been swapped according to WHICH.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.
*> \endverbatim
      SUBROUTINE SLAREV( WHICH, M, N, A, LDA )
*     .. Scalar Arguments ..
      CHARACTER          WHICH
      INTEGER            M, N, LDA
*     .. Array Arguments ..
      REAL   A( LDA, * )
*     .. Local Scalars ..
      INTEGER            I, J, M2, N2
      LOGICAL            ROWODD, COLODD
      REAL   TMPA, TMPB, TMPC, TMPD
*     .. Executable Statements ..

      M2 = M / 2
      N2 = N / 2
      ROWODD = MOD(M, 2) .NE. 0
      COLODD = MOD(N, 2) .NE. 0

      IF ( WHICH .EQ. 'B' ) THEN
         DO I = 1, N2
            DO J = 1, M2
               TMPA = A( J, I )
               TMPB = A( J, N - I + 1 )
               TMPC = A( M - J + 1, I )
               TMPD = A( M - J + 1, N - I + 1 )
               
               A( J, I ) = TMPD
               A( J, N - I + 1 ) = TMPC
               A( M - J + 1, I ) = TMPB
               A( M - J + 1, N - I + 1 ) = TMPA
            END DO
            IF ( ROWODD ) THEN
               TMPA = A( M2 + 1, I )
               A( M2 + 1, I ) = A( M2 + 1, N - I + 1 )
               A( M2 + 1, N - I + 1 ) = TMPA
            END IF
         END DO
         IF ( COLODD ) THEN
            DO J = 1, M2
               TMPA = A( J, N2 + 1 )
               A( J, N2 + 1 ) = A( M - J + 1, N2 + 1 )
               A( M - J + 1, N2 + 1 ) = TMPA
            END DO
         END IF

      ELSE IF ( WHICH .EQ. 'R' ) THEN
         DO I = 1, N
            DO J = 1, M2
               TMPA = A( J, I )
               A( J, I ) = A( M - J + 1, I )
               A( M - J + 1, I ) = TMPA
            END DO
         END DO

      ELSE IF ( WHICH .EQ. 'C' ) THEN
         DO I = 1, N2
            DO J = 1, M
               TMPA = A( J, I )
               A( J, I ) = A( J, N - I + 1 )
               A( J, N - I + 1 ) = TMPA
            END DO
         END DO
      END IF

      RETURN
      END