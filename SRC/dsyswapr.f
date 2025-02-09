*> \brief \b DSYSWAPR applies an elementary permutation on the rows and columns of a symmetric matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download DSYSWAPR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dsyswapr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dsyswapr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dsyswapr.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSYSWAPR( UPLO, N, A, LDA, I1, I2)
*
*       .. Scalar Arguments ..
*       CHARACTER        UPLO
*       INTEGER          I1, I2, LDA, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A( LDA, N )
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DSYSWAPR applies an elementary permutation on the rows and the columns of
*> a symmetric matrix.
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
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix A.  N >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,*)
*>          On entry, the N-by-N matrix A. On exit, the permuted matrix
*>          where the rows I1 and I2 and columns I1 and I2 are interchanged.
*>          If UPLO = 'U', the interchanges are applied to the upper
*>          triangular part and the strictly lower triangular part of A is
*>          not referenced; if UPLO = 'L', the interchanges are applied to
*>          the lower triangular part and the part of A above the diagonal
*>          is not referenced.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,N).
*> \endverbatim
*>
*> \param[in] I1
*> \verbatim
*>          I1 is INTEGER
*>          Index of the first row to swap
*> \endverbatim
*>
*> \param[in] I2
*> \verbatim
*>          I2 is INTEGER
*>          Index of the second row to swap
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
*> \ingroup heswapr
*
*  =====================================================================
      SUBROUTINE DSYSWAPR( UPLO, N, A, LDA, I1, I2)
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER        UPLO
      INTEGER          I1, I2, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A( LDA, * )
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      DOUBLE PRECISION   TMP
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL         DSWAP
*     ..
*     .. Executable Statements ..
*
      UPPER = LSAME( UPLO, 'U' )
      IF (UPPER) THEN
*
*         UPPER
*         first swap
*          - swap column I1 and I2 from I1 to I1-1
         CALL DSWAP( I1-1, A(1,I1), 1, A(1,I2), 1 )
*
*          second swap :
*          - swap A(I1,I1) and A(I2,I2)
*          - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
         TMP=A(I1,I1)
         A(I1,I1)=A(I2,I2)
         A(I2,I2)=TMP
*
         CALL DSWAP( I2-I1-1, A(I1,I1+1), LDA, A(I1+1,I2), 1 )
*
*          third swap
*          - swap row I1 and I2 from I2+1 to N
         IF ( I2.LT.N )
     $      CALL DSWAP( N-I2, A(I1,I2+1), LDA, A(I2,I2+1), LDA )
*
        ELSE
*
*         LOWER
*         first swap
*          - swap row I1 and I2 from I1 to I1-1
         CALL DSWAP( I1-1, A(I1,1), LDA, A(I2,1), LDA )
*
*         second swap :
*          - swap A(I1,I1) and A(I2,I2)
*          - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
          TMP=A(I1,I1)
          A(I1,I1)=A(I2,I2)
          A(I2,I2)=TMP
*
          CALL DSWAP( I2-I1-1, A(I1+1,I1), 1, A(I2,I1+1), LDA )
*
*         third swap
*          - swap col I1 and I2 from I2+1 to N
         IF ( I2.LT.N )
     $      CALL DSWAP( N-I2, A(I2+1,I1), 1, A(I2+1,I2), 1 )
*
      ENDIF
      END SUBROUTINE DSYSWAPR

