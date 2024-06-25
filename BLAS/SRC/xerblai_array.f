*> \brief \b XERBLAI_ARRAY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE XERBLAI_ARRAY(SRNAME_ARRAY, SRNAME_LEN, INFO, INDX)
*
*       .. Scalar Arguments ..
*       INTEGER SRNAME_LEN, INFO, INDX
*       ..
*       .. Array Arguments ..
*       CHARACTER(1) SRNAME_ARRAY(SRNAME_LEN)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> XERBLAI_ARRAY assists other languages in calling XERBLAI, the LAPACK
*> and BLAS error handler.  Rather than taking a Fortran string argument
*> as the function's name, XERBLAI_ARRAY takes an array of single
*> characters along with the array's length.  XERBLAI_ARRAY then copies
*> up to 32 characters of that array into a Fortran string and passes
*> that to XERBLAI.  If called with a non-positive SRNAME_LEN,
*> XERBLAI_ARRAY will call XERBLAI with a string of all blank characters.
*>
*> Say some macro or other device makes XERBLAI_ARRAY available to C99
*> by a name lapack_xerbla and with a common Fortran calling convention.
*> Then a C99 program could invoke XERBLAI via:
*>    {
*>      int flen = strlen(__func__);
*>      lapack_xerblai(__func__, &flen, &info, &indx);
*>    }
*>
*> Providing XERBLAI_ARRAY is not necessary for intercepting LAPACK
*> errors.  XERBLAI_ARRAY calls XERBLAI.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SRNAME_ARRAY
*> \verbatim
*>          SRNAME_ARRAY is CHARACTER(1) array, dimension (SRNAME_LEN)
*>          The name of the routine which called XERBLAI_ARRAY.
*> \endverbatim
*>
*> \param[in] SRNAME_LEN
*> \verbatim
*>          SRNAME_LEN is INTEGER
*>          The length of the name in SRNAME_ARRAY.
*> \endverbatim
*>
*> \param[in] INFO
*> \verbatim
*>          INFO is INTEGER
*>          The position of the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*>
*> \param[in] INDX
*> \verbatim
*>          INDX is INTEGER
*>          The position at the invalid parameter in the parameter list
*>          of the calling routine.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*> \author Igor S. Gerasimov
*
*> \ingroup xerbla_array
*
*  =====================================================================
      SUBROUTINE XERBLAI_ARRAY(SRNAME_ARRAY, SRNAME_LEN, INFO, INDX)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER SRNAME_LEN, INFO, INDX
*     ..
*     .. Array Arguments ..
      CHARACTER(1) SRNAME_ARRAY(SRNAME_LEN)
*     ..
*
* =====================================================================
*
*     ..
*     .. Local Scalars ..
      INTEGER I
*     ..
*     .. Local Arrays ..
      CHARACTER*32 SRNAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MIN, LEN
*     ..
*     .. External Functions ..
      EXTERNAL XERBLAI
*     ..
*     .. Executable Statements ..
      SRNAME = ' '
      DO I = 1, MIN( SRNAME_LEN, LEN( SRNAME ) )
         SRNAME( I:I ) = SRNAME_ARRAY( I )
      END DO

      CALL XERBLAI( SRNAME, INFO, INDX )

      RETURN
*
*     End of XERBLAI_ARRAY
*
      END
