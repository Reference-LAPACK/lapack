*> \brief \b CLAEV2 computes the eigenvalues and eigenvectors of a 2-by-2 symmetric/Hermitian matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> Download CLAEV2 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claev2.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claev2.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claev2.f">
*> [TXT]</a>
*
*  Definition:
*  ===========
*
*       SUBROUTINE CLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*       .. Scalar Arguments ..
*       REAL               CS1, RT1, RT2
*       COMPLEX            A, B, C, SN1
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CLAEV2 computes the eigendecomposition of a 2-by-2 Hermitian matrix
*>    [  A         B  ]
*>    [  CONJG(B)  C  ].
*> On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*> eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*> eigenvector for RT1, giving the decomposition
*>
*> [ CS1  CONJG(SN1) ] [    A     B ] [ CS1 -CONJG(SN1) ] = [ RT1  0  ]
*> [-SN1     CS1     ] [ CONJG(B) C ] [ SN1     CS1     ]   [  0  RT2 ].
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] A
*> \verbatim
*>          A is COMPLEX
*>         The (1,1) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] B
*> \verbatim
*>          B is COMPLEX
*>         The (1,2) element and the conjugate of the (2,1) element of
*>         the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is COMPLEX
*>         The (2,2) element of the 2-by-2 matrix.
*> \endverbatim
*>
*> \param[out] RT1
*> \verbatim
*>          RT1 is REAL
*>         The eigenvalue of larger absolute value.
*> \endverbatim
*>
*> \param[out] RT2
*> \verbatim
*>          RT2 is REAL
*>         The eigenvalue of smaller absolute value.
*> \endverbatim
*>
*> \param[out] CS1
*> \verbatim
*>          CS1 is REAL
*> \endverbatim
*>
*> \param[out] SN1
*> \verbatim
*>          SN1 is COMPLEX
*>         The vector (CS1, SN1) is a unit right eigenvector for RT1.
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
*> \ingroup laev2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  RT1 is accurate to a few ulps barring over/underflow.
*>
*>  RT2 may be inaccurate if there is massive cancellation in the
*>  determinant A*C-B*B; higher precision or correctly rounded or
*>  correctly truncated arithmetic would be needed to compute RT2
*>  accurately in all cases.
*>
*>  CS1 and SN1 are accurate to a few ulps barring over/underflow.
*>
*>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*>  Underflow is harmless if the input data is 0 or exceeds
*>     underflow_threshold / macheps.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CLAEV2( A, B, C, RT1, RT2, CS1, SN1 )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL               CS1, RT1, RT2
      COMPLEX            A, B, C, SN1
*     ..
*
* =====================================================================
*
*     .. Parameters ..
      REAL               ZERO
      PARAMETER          ( ZERO = 0.0E0 )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E0 )
*     ..
*     .. Local Scalars ..
      REAL               T
      COMPLEX            W
*     ..
*     .. External Subroutines ..
      EXTERNAL           SLAEV2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, CONJG, REAL
*     ..
*     .. Executable Statements ..
*
      IF( ABS( B ).EQ.ZERO ) THEN
         W = ONE
      ELSE
         W = CONJG( B ) / ABS( B )
      END IF
      CALL SLAEV2( REAL( A ), ABS( B ), REAL( C ), RT1, RT2, CS1, T )
      SN1 = W*T
      RETURN
*
*     End of CLAEV2
*
      END
