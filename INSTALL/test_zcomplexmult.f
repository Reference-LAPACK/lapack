*> \brief zmul tests the robustness and precision of the double complex multiplication
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Authors:
*  ========
*
*> \author Weslley S. Pereira, University of Colorado Denver, U.S.
*
*> \verbatim
*>
*> Tests:
*>
*> (a) Inf inputs:
*>    (1) y = ( Inf + 0   * I)
*>    (2) y = (-Inf + 0   * I)
*>    (3) y = ( 0   + Inf * I)
*>    (4) y = ( 0   - Inf * I)
*>    (5) y = ( Inf + Inf * I)
*> Tests:
*>    (a) 0 * y is NaN.
*>    (b) 1 * y is y is either y or NaN.
*>    (c) y * y is either  Inf or NaN (cases 1 and 3),
*>                 either -Inf or NaN (cases 2 and 4),
*>                 NaN (case 5).
*>
*> (b) NaN inputs:
*>    (1) y = (NaN + 0   * I)
*>    (2) y = (0   + NaN * I)
*>    (3) y = (NaN + NaN * I)
*> Tests:
*>    (a) 0 * y is NaN.
*>    (b) 1 * y is NaN.
*>    (c) y * y is NaN.
*>
*> \endverbatim
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      program zmul
*
*  -- LAPACK test routine --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

*     ..
*     .. Constants ..
      integer           DP
      parameter       ( DP = KIND( 0.0d0 ) )
      integer           nNaN, nInf
      parameter       ( nNaN = 3, nInf = 5 )
      double complex    czero, cone
      parameter       ( czero = CMPLX( 0.0d0, 0.0d0, KIND=DP ),
     $                  cone  = CMPLX( 1.0d0, 0.0d0, KIND=DP ) )
*     ..
*     .. Local Variables ..
      integer           i, nFailingTests, nTests
      double precision  aInf, aNaN, OV
      double complex    Y, R, cInf( nInf ), cNaN( nNaN )
*
*     .. Intrinsic Functions ..
      intrinsic         CMPLX, HUGE

*
*     .. Initialize error counts ..
      nFailingTests = 0
      nTests = 0
*
*     .. Inf entries ..
      OV = HUGE(0.0d0)
      aInf = OV * 2
      cInf(1) = CMPLX( aInf, 0.0d0, KIND=DP )
      cInf(2) = CMPLX(-aInf, 0.0d0, KIND=DP )
      cInf(3) = CMPLX( 0.0d0, aInf, KIND=DP )
      cInf(4) = CMPLX( 0.0d0,-aInf, KIND=DP )
      cInf(5) = CMPLX( aInf,  aInf, KIND=DP )
*
*     .. NaN entries ..
      aNaN = aInf / aInf
      cNaN(1) = CMPLX( aNaN, 0.0d0, KIND=DP )
      cNaN(2) = CMPLX( 0.0d0, aNaN, KIND=DP )
      cNaN(3) = CMPLX( aNaN,  aNaN, KIND=DP )

*
*     .. Tests ..
*
*     Test (a) Infs
      do 10 i = 1, nInf
          nTests = nTests + 3
          Y = cInf(i)
          R = czero * Y
          if( R .eq. R ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'ia',i, czero, Y, R, 'NaN'
          endif
          R = cone * Y
          if( (R .ne. Y) .and. (R .eq. R) ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'ib',i, cone, Y, R,
     $                               'the input and NaN'
          endif
          R = Y * Y
          if( (i.eq.1) .or. (i.eq.2) ) then
              if( (R .ne. cInf(1)) .and. (R .eq. R) ) then
                  nFailingTests = nFailingTests + 1
                  WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, 'Inf and NaN'
              endif
          else if( (i.eq.3) .or. (i.eq.4) ) then
              if( (R .ne. cInf(2)) .and. (R .eq. R) ) then
                  nFailingTests = nFailingTests + 1
                  WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, '-Inf and NaN'
              endif
          else 
              if( R .eq. R ) then
                  nFailingTests = nFailingTests + 1
                  WRITE( *, FMT = 9998 ) 'ic',i, Y, Y, R, 'NaN'
              endif
          endif
  10  continue
*
*     Test (b) NaNs
      do 20 i = 1, nNaN
          nTests = nTests + 3
          Y = cNaN(i)
          R = czero * Y
          if( R .eq. R ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'na',i, czero, Y, R, 'NaN'
          endif
          R = cone * Y
          if( R .eq. R ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'nb',i, cone, Y, R, 'NaN'
          endif
          R = Y * Y
          if( R .eq. R ) then
              nFailingTests = nFailingTests + 1
              WRITE( *, FMT = 9998 ) 'nc',i, Y, Y, R, 'NaN'
          endif
  20  continue
*
      if( nFailingTests .gt. 0 ) then
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests,
     $      " pass for complex multiplication,", nFailingTests," fail."
      else
         print *, "# All tests pass for complex multiplication."
      endif
*
*     .. Formats ..
 9998 FORMAT( '[',A2,I1, '] (', (ES24.16E3,SP,ES24.16E3,"*I"), ') * (',
     $         (ES24.16E3,SP,ES24.16E3,"*I"), ') = (',
     $         (ES24.16E3,SP,ES24.16E3,"*I"), ') differs from ', A17 )
*
*     End of zmul
*
      END