*> \brief \b SKYCONV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download SKYCONV + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/skyconv.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/skyconv.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/skyconv.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE SKYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*       .. Scalar Arguments ..
*       CHARACTER          UPLO, WAY
*       INTEGER            INFO, LDA, N
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       REAL               A( LDA, * ), E( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> SKYCONV convert A given by TRF into L and D and vice-versa.
*> Get Non-diag elements of D (returned in workspace) and
*> apply or reverse permutation done in TRF.
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
*>          = 'U':  Upper triangular, form is A = U*D*U**T;
*>          = 'L':  Lower triangular, form is A = L*D*L**T.
*> \endverbatim
*>
*> \param[in] WAY
*> \verbatim
*>          WAY is CHARACTER*1
*>          = 'C': Convert
*>          = 'R': Revert
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is REAL array, dimension (LDA,N)
*>          The block diagonal matrix D and the multipliers used to
*>          obtain the factor U or L as computed by SKYTRF.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] IPIV
*> \verbatim
*>          IPIV is INTEGER array, dimension (N)
*>          Details of the interchanges and the block structure of D
*>          as determined by SKYTRF.
*> \endverbatim
*>
*> \param[out] E
*> \verbatim
*>          E is REAL array, dimension (N)
*>          E stores the supdiagonal/subdiagonal of the skew-symmetric
*>          2-by-2 block diagonal matrix D in LDLT.
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
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup kyconv
*
*  =====================================================================
      SUBROUTINE SKYCONV( UPLO, WAY, N, A, LDA, IPIV, E, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, WAY
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      REAL               A( LDA, * ), E( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E+0 )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Local Scalars ..
      LOGICAL            UPPER, CONVERT
      INTEGER            I, IP, J
      REAL               TEMP
*     ..
*     .. Executable Statements ..
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      CONVERT = LSAME( WAY, 'C' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.CONVERT .AND. .NOT.LSAME( WAY, 'R' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5

      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SKYCONV', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*      A is UPPER
*
*      Convert A (A is upper)
*
*        Convert VALUE
*
         IF ( CONVERT ) THEN
            I=N
            E(1)=ZERO
            DO WHILE ( I .GT. 1 )
               E(I)=A(I-1,I)
               A(I-1,I)=ZERO
               I=I-2
            END DO
*
*        Convert PERMUTATIONS
*
         I=N-2
         DO WHILE ( I .GT. 1 )
            IF( IPIV(I) .GT. 0) THEN
               IP=IPIV(I)
                  DO 12 J= I+1,N
                     TEMP=A(IP,J)
                     A(IP,J)=A(I-1,J)
                     A(I-1,J)=TEMP
 12            CONTINUE
            ELSEIF( IPIV(I) .LT. 0) THEN
               IP=-IPIV(I)
				   DO 13 J= I+1,N
                  TEMP=A(I,J)
                  A(I,J)=A(I-1,J)
                  A(I-1,J)=TEMP
                  
                  TEMP=A(IP,J)
                  A(IP,J)=A(I-1,J)
                  A(I-1,J)=TEMP
 13            CONTINUE
            ENDIF
            I=I-2
         END DO

         ELSE
*
*      Revert A (A is upper)
*
*
*        Revert PERMUTATIONS
*
            I=2
            DO WHILE ( I .LT. N-1 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  DO J= I+1,N
                     TEMP=A(IP,J)
                     A(IP,J)=A(I-1,J)
                     A(I-1,J)=TEMP
                  END DO
               ELSEIF( IPIV(I) .LT. 0 ) THEN
                  IP=-IPIV(I)
                  DO J= I+1,N
                     TEMP=A(IP,J)
                     A(IP,J)=A(I-1,J)
                     A(I-1,J)=TEMP

                     TEMP=A(I,J)
                     A(I,J)=A(I-1,J)
                     A(I-1,J)=TEMP
                 END DO
               ENDIF
               I=I+2
            END DO
*
*        Revert VALUE
*
            I=N
            DO WHILE ( I .GT. 1 )
               A(I-1,I)=E(I)
               I=I-2
            END DO
         END IF
      ELSE
*
*      A is LOWER
*
         IF ( CONVERT ) THEN
*
*      Convert A (A is lower)
*
*
*        Convert VALUE
*
            I=1
            E(N)=ZERO
            DO WHILE ( I .LT. N )
               E(I)=A(I+1,I)
               A(I+1,I)=ZERO
               I=I+2
            END DO
*
*        Convert PERMUTATIONS
*
         I=3
         DO WHILE ( I .LT. N )
            IF( IPIV(I) .GT. 0 ) THEN
               IP=IPIV(I)
               DO 22 J= 1,I-1
                  TEMP=A(IP,J)
                  A(IP,J)=A(I+1,J)
                  A(I+1,J)=TEMP
 22            CONTINUE
            ELSEIF( IPIV(I) .LT. 0 ) THEN
               IP=-IPIV(I)
               DO 23 J= 1,I-1
                  TEMP=A(I,J)
                  A(I,J)=A(I+1,J)
                  A(I+1,J)=TEMP

                  TEMP=A(IP,J)
                  A(IP,J)=A(I+1,J)
                  A(I+1,J)=TEMP
 23           CONTINUE
           ENDIF
           I=I+2
         END DO
         ELSE
*
*      Revert A (A is lower)
*
*
*        Revert PERMUTATIONS
*
            I=N-1
            DO WHILE ( I .GT. 2 )
               IF( IPIV(I) .GT. 0 ) THEN
                  IP=IPIV(I)
                  DO J= 1,I-1
                     TEMP=A(I+1,J)
                     A(I+1,J)=A(IP,J)
                     A(IP,J)=TEMP
                  END DO
               ELSEIF( IPIV(I) .LT. 0 ) THEN
                  IP=-IPIV(I)
                  DO J= 1,I-1
                     TEMP=A(I+1,J)
                     A(I+1,J)=A(IP,J)
                     A(IP,J)=TEMP

                     TEMP=A(I+1,J)
                     A(I+1,J)=A(I,J)
                     A(I,J)=TEMP
                  END DO
               ENDIF
               I=I-2
            END DO
*
*        Revert VALUE
*
            I=1
            DO WHILE ( I .LT. N )
               A(I+1,I)=E(I)
               I=I+2
            END DO
         END IF
      END IF

      RETURN
*
*     End of SKYCONV
*
      END
