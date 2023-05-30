*> \brief zabs tests the robustness and precision of the intrinsic ABS for double complex
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \author Weslley S. Pereira, University of Colorado Denver, U.S.
*
*> \verbatim
*>
*> Real values for test:
*> (1) x = 2**m, where m = MINEXPONENT-DIGITS, ..., MINEXPONENT-1. Stop on the first success.
*>     Mind that not all platforms might implement subnormal numbers.
*> (2) x = 2**m, where m = MINEXPONENT, ..., 0. Stop on the first success.
*> (3) x = OV, where OV is the overflow threshold. OV^2 overflows but the norm is OV.
*> (4) x = 2**m, where m = MAXEXPONENT-1, ..., 1. Stop on the first success.
*>
*> Tests:
*> (a) y = x + 0 * I, |y| = x
*> (b) y = 0 + x * I, |y| = x
*> (c) y = (3/4)*x + x * I, |y| = (5/4)*x whenever (3/4)*x and (5/4)*x can be exactly stored
*> (d) y = (1/2)*x + (1/2)*x * I, |y| = (1/2)*x*sqrt(2) whenever (1/2)*x can be exactly stored
*>
*> Special cases:
*>
*> (i) Inf propagation
*>    (1) y = Inf + 0 * I, |y| is Inf.
*>    (2) y =-Inf + 0 * I, |y| is Inf.
*>    (3) y = 0 + Inf * I, |y| is Inf.
*>    (4) y = 0 - Inf * I, |y| is Inf.
*>    (5) y = Inf + Inf * I, |y| is Inf.
*>
*> (n) NaN propagation
*>    (1) y = NaN + 0 * I, |y| is NaN.
*>    (2) y = 0 + NaN * I, |y| is NaN.
*>    (3) y = NaN + NaN * I, |y| is NaN.
*>
*> \endverbatim
*
*> \ingroup auxOTHERauxiliary
*
*  =====================================================================
      program zabs
*
*  -- LAPACK test routine --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..

*     ..
*     .. Local parameters ..
      logical           debug
      parameter       ( debug = .false. )
      integer           N, nNaN, nInf
      parameter       ( N = 4, nNaN = 3, nInf = 5 )
      double precision  threeFourth, fiveFourth, oneHalf
      parameter       ( threeFourth = 3.0d0 / 4,
     $                  fiveFourth = 5.0d0 / 4,
     $                  oneHalf = 1.0d0 / 2 )
*     ..
*     .. Local Variables ..
      integer           i, min, Max, m, subnormalTreatedAs0,
     $                  caseAFails, caseBFails, caseCFails, caseDFails,
     $                  caseEFails, caseFFails, nFailingTests, nTests
      double precision  X( N ), R, answerC,
     $                  answerD, aInf, aNaN, relDiff, b,
     $                  eps, blueMin, blueMax, Xj, stepX(N), limX(N)
      double complex    Y, cInf( nInf ), cNaN( nNaN )
*
*     .. Intrinsic Functions ..
      intrinsic         ABS, DBLE, RADIX, CEILING, TINY, DIGITS, SQRT,
     $                  MAXEXPONENT, MINEXPONENT, FLOOR, HUGE, DCMPLX,
     $                  EPSILON

*
*     .. Initialize error counts ..
      subnormalTreatedAs0 = 0
      caseAFails = 0
      caseBFails = 0
      caseCFails = 0
      caseDFails = 0
      caseEFails = 0
      caseFFails = 0
      nFailingTests = 0
      nTests = 0
*
*     .. Initialize machine constants ..
      min = MINEXPONENT(0.0d0)
      Max = MAXEXPONENT(0.0d0)
      m = DIGITS(0.0d0)
      b = DBLE(RADIX(0.0d0))
      eps = EPSILON(0.0d0)
      blueMin = b**CEILING( (min - 1) * 0.5d0 )
      blueMax = b**FLOOR( (Max - m + 1) * 0.5d0 )
*
*     .. Vector X ..
      X(1) = TINY(0.0d0) * b**( DBLE(1-m) )
      X(2) = TINY(0.0d0)
      X(3) = HUGE(0.0d0)
      X(4) = b**( DBLE(Max-1) )
*
*     .. Then modify X using the step ..
      stepX(1) = 2.0
      stepX(2) = 2.0
      stepX(3) = 0.0
      stepX(4) = 0.5
*
*     .. Up to the value ..
      limX(1) = X(2)
      limX(2) = 1.0
      limX(3) = 0.0
      limX(4) = 2.0
*
*     .. Inf entries ..
      aInf = X(3) * 2
      cInf(1) = DCMPLX( aInf, 0.0d0 )
      cInf(2) = DCMPLX(-aInf, 0.0d0 )
      cInf(3) = DCMPLX( 0.0d0, aInf )
      cInf(4) = DCMPLX( 0.0d0,-aInf )
      cInf(5) = DCMPLX( aInf,  aInf )
*
*     .. NaN entries ..
      aNaN = aInf / aInf
      cNaN(1) = DCMPLX( aNaN, 0.0d0 )
      cNaN(2) = DCMPLX( 0.0d0, aNaN )
      cNaN(3) = DCMPLX( aNaN,  aNaN )

*
*     .. Tests ..
*
      if( debug ) then
        print *, '# X :=', X
        print *, '# Blue min constant :=', blueMin
        print *, '# Blue max constant :=', blueMax
      endif
*
      Xj = X(1)
      if( Xj .eq. 0.0d0 ) then
        subnormalTreatedAs0 = subnormalTreatedAs0 + 1
        if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
            print *, "!! fl( subnormal ) may be 0"
        endif
      else
        do 100 i = 1, N
            Xj = X(i)
            if( Xj .eq. 0.0d0 ) then
                subnormalTreatedAs0 = subnormalTreatedAs0 + 1
                if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
                    print *, "!! fl( subnormal ) may be 0"
                endif
            endif
 100    continue
      endif
*
*     Test (a) y = x + 0 * I, |y| = x
      do 10 i = 1, N
        Xj = X(i)
        if( Xj .eq. 0.0d0 ) then
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
                print *, "!! [a] fl( subnormal ) may be 0"
            endif
        else
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( Xj, 0.0d0 )
                R = ABS( Y )
                if( R .ne. Xj ) then
                    caseAFails = caseAFails + 1
                    if( caseAFails .eq. 1 ) then
                        print *, "!! Some ABS(x+0*I) differ from ABS(x)"
                    endif
                    WRITE( 0, FMT = 9999 ) 'a',i, Xj, '(1+0*I)', R, Xj
                endif
                Xj = Xj * stepX(i)
            end do
        endif
  10  continue
*
*     Test (b) y = 0 + x * I, |y| = x
      do 20 i = 1, N
        Xj = X(i)
        if( Xj .eq. 0.0d0 ) then
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
                print *, "!! [b] fl( subnormal ) may be 0"
            endif
        else
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                Y = DCMPLX( 0.0d0, Xj )
                R = ABS( Y )
                if( R .ne. Xj ) then
                    caseBFails = caseBFails + 1
                    if( caseBFails .eq. 1 ) then
                        print *, "!! Some ABS(0+x*I) differ from ABS(x)"
                    endif
                    WRITE( 0, FMT = 9999 ) 'b',i, Xj, '(0+1*I)', R, Xj
                endif
                Xj = Xj * stepX(i)
            end do
        endif
  20  continue
*
*     Test (c) y = (3/4)*x + x * I, |y| = (5/4)*x
      do 30 i = 1, N
        if( i .eq. 3 ) go to 30
        if( i .eq. 1 ) then
            Xj = 4*X(i)
        else
            Xj = X(i)
        endif
        if( Xj .eq. 0.0d0 ) then
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
                print *, "!! [c] fl( subnormal ) may be 0"
            endif
        else
            do while( Xj .ne. limX(i) )
                nTests = nTests + 1
                answerC = fiveFourth * Xj
                Y = DCMPLX( threeFourth * Xj, Xj )
                R = ABS( Y )
                if( R .ne. answerC ) then
                    caseCFails = caseCFails + 1
                    if( caseCFails .eq. 1 ) then
                        print *, 
     $              "!! Some ABS(x*(3/4+I)) differ from (5/4)*ABS(x)"
                    endif
                    WRITE( 0, FMT = 9999 ) 'c',i, Xj, '(3/4+I)', R,
     $                                      answerC
                endif
                Xj = Xj * stepX(i)
            end do
        endif
  30  continue
*
*     Test (d) y = (1/2)*x + (1/2)*x * I, |y| = (1/2)*x*sqrt(2)
      do 40 i = 1, N
        if( i .eq. 1 ) then
            Xj = 2*X(i)
        else
            Xj = X(i)
        endif
        if( Xj .eq. 0.0d0 ) then
            subnormalTreatedAs0 = subnormalTreatedAs0 + 1
            if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
                print *, "!! [d] fl( subnormal ) may be 0"
            endif
        else
            do while( Xj .ne. limX(i) )
                answerD = (oneHalf * Xj) * SQRT(2.0d0)
                if( answerD .eq. 0.0d0 ) then
                    subnormalTreatedAs0 = subnormalTreatedAs0 + 1
                    if( debug .or. subnormalTreatedAs0 .eq. 1 ) then
                        print *, "!! [d] fl( subnormal ) may be 0"
                    endif
                else
                    nTests = nTests + 1
                    Y = DCMPLX( oneHalf * Xj, oneHalf * Xj )
                    R = ABS( Y )
                    relDiff = ABS(R-answerD)/answerD
                    if( relDiff .ge. (0.5*eps) ) then
                        caseDFails = caseDFails + 1
                        if( caseDFails .eq. 1 ) then
                            print *, 
     $                "!! Some ABS(x*(1+I)) differ from sqrt(2)*ABS(x)"
                        endif
                        WRITE( 0, FMT = 9999 ) 'd',i, (oneHalf*Xj),
     $                                  '(1+1*I)', R, answerD
                    endif
                endif
                Xj = Xj * stepX(i)
            end do
        endif 
  40  continue
*
*     Test (e) Infs
      do 50 i = 1, nInf
        nTests = nTests + 1
        Y = cInf(i)
        R = ABS( Y )
        if( .not.(R .gt. HUGE(0.0d0)) ) then
            caseEFails = caseEFails + 1
            WRITE( *, FMT = 9997 ) 'i',i, Y, R
        endif
  50  continue
*
*     Test (f) NaNs
      do 60 i = 1, nNaN
        nTests = nTests + 1
        Y = cNaN(i)
        R = ABS( Y )
        if( R .eq. R ) then
            caseFFails = caseFFails + 1
            WRITE( *, FMT = 9998 ) 'n',i, Y, R
        endif
  60  continue
*
*     If any test fails, displays a message
      nFailingTests = caseAFails + caseBFails + caseCFails + caseDFails
     $                + caseEFails + caseFFails
      if( nFailingTests .gt. 0 ) then
         print *, "# ", nTests-nFailingTests, " tests out of ", nTests,
     $            " pass for ABS(a+b*I),", nFailingTests, " tests fail."
      else
         print *, "# All tests pass for ABS(a+b*I)"
      endif
*
*     If anything was written to stderr, print the message
      if( (caseAFails .gt. 0) .or. (caseBFails .gt. 0) .or.
     $    (caseCFails .gt. 0) .or. (caseDFails .gt. 0) ) then
         print *, "# Please check the failed ABS(a+b*I) in [stderr]"
      endif
*
*     .. Formats ..
 9997 FORMAT( '[',A1,I1, '] ABS(', (ES8.1,SP,ES8.1,"*I"), ' ) = ',
     $        ES8.1, ' differs from Inf' )
*
 9998 FORMAT( '[',A1,I1, '] ABS(', (ES8.1,SP,ES8.1,"*I"), ' ) = ',
     $        ES8.1, ' differs from NaN' )
*
 9999 FORMAT( '[',A1,I1, '] ABS(', ES24.16E3, ' * ', A7, ' ) = ',
     $         ES24.16E3, ' differs from ', ES24.16E3 )
*
*     End of zabs
*
      END
