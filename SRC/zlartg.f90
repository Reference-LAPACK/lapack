!> \brief \b ZLARTG generates a plane rotation with real cosine and complex sine.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARTG( F, G, C, S, R )
!
!       .. Scalar Arguments ..
!       REAL(wp)           C
!       COMPLEX(wp)        F, G, R, S
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARTG generates a plane rotation so that
!>
!>    [  C         S  ] . [ F ]  =  [ R ]
!>    [ -conjg(S)  C  ]   [ G ]     [ 0 ]
!>
!> where C is real and C**2 + |S|**2 = 1.
!>
!> The mathematical formulas used for C and S are
!>
!>    sgn(x) = {  x / |x|,   x != 0
!>             {  1,         x = 0
!>
!>    R = sgn(F) * sqrt(|F|**2 + |G|**2)
!>
!>    C = |F| / sqrt(|F|**2 + |G|**2)
!>
!>    S = sgn(F) * conjg(G) / sqrt(|F|**2 + |G|**2)
!>
!> When F and G are real, the formulas simplify to C = F/R and
!> S = G/R, and the returned values of C, S, and R should be
!> identical to those returned by DLARTG.
!>
!> The algorithm used to compute these quantities incorporates scaling
!> to avoid overflow or underflow in computing the square root of the
!> sum of squares.
!>
!> This is a faster version of the BLAS1 routine ZROTG, except for
!> the following differences:
!>    F and G are unchanged on return.
!>    If G=0, then C=1 and S=0.
!>    If F=0, then C=0 and S is chosen so that R is real.
!>
!> Below, wp=>dp stands for double precision from LA_CONSTANTS module.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is COMPLEX(wp)
!>          The first component of vector to be rotated.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is COMPLEX(wp)
!>          The second component of vector to be rotated.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL(wp)
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX(wp)
!>          The sine of the rotation.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX(wp)
!>          The nonzero component of the rotated vector.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Weslley Pereira, University of Colorado Denver, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Anderson E. (2017)
!>  Algorithm 978: Safe Scaling in the Level 1 BLAS
!>  ACM Trans Math Softw 44:1--28
!>  https://doi.org/10.1145/3061665
!>
!> \endverbatim
!
subroutine ZLARTG( f, g, c, s, r )
   use LA_CONSTANTS, &
   only: wp=>dp, zero=>dzero, one=>done, two=>dtwo, czero=>zzero, &
         rtmin=>drtmin, rtmax=>drtmax, safmin=>dsafmin, safmax=>dsafmax
!
!  -- LAPACK auxiliary routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2021
!
!  .. Scalar Arguments ..
   real(wp)           c
   complex(wp)        f, g, r, s
!  ..
!  .. Local Scalars ..
   real(wp) :: d, f1, f2, g1, g2, h2, w2, u, uu, v, vv, w
   complex(wp) :: fs, gs, t
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, conjg, max, min, real, sqrt
!  ..
!  .. Statement Functions ..
   real(wp) :: ABSSQ
!  ..
!  .. Statement Function definitions ..
   ABSSQ( t ) = real( t )**2 + aimag( t )**2
!  ..
!  .. Executable Statements ..
!
   if( g == czero ) then
      c = one
      s = czero
      r = f
   else if( f == czero ) then
      c = zero
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         d = abs( g )
         s = conjg( g ) / d
         r = d
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, g1 ) )
         uu = one / u
         gs = g*uu
         d = abs( g2 )
         s = conjg( gs ) / d
         r = d*u
      end if
   else
      f1 = max( abs(real(f)), abs(aimag(f)) )
      g1 = max( abs(real(g)), abs(aimag(g)) )
      if( f1 > rtmin .and. f1 < rtmax .and. &
          g1 > rtmin .and. g1 < rtmax ) then
!
!        Use unscaled algorithm
!
         f2 = ABSSQ( f )
         g2 = ABSSQ( g )
         h2 = f2 + g2
         d = sqrt( one + ( g2/f2 ) )
         r = f*d
         c = one / d
         s = conjg( g )*( r / h2 )
      else
!
!        Use scaled algorithm
!
         u = min( safmax, max( safmin, f1, g1 ) )
         uu = one / u
         gs = g*uu
         g2 = ABSSQ( gs )
         if( f1*uu < rtmin ) then
!
!           f is not well-scaled when scaled by g1.
!           Use a different scaling for f.
!
            v = min( safmax, max( safmin, f1 ) )
            vv = one / v
            w = v * uu
            w2 = w**2
            fs = f*vv
            f2 = ABSSQ( fs )
            h2 = f2*w2 + g2
         else
!
!           Otherwise use the same scaling for f and g.
!
            w = one
            w2 = one
            fs = f*uu
            f2 = ABSSQ( fs )
            h2 = f2 + g2
         end if
         d = sqrt( w2 + ( g2/f2 ) )
         c = w / d
         r = fs*d
         s = conjg( gs )*( r / h2 )
         r = r*u
      end if
   end if
   return
end subroutine
