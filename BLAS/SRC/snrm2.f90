!> \brief \b SNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    SNRM2 := sqrt( x'*x ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER, storage spacing between elements of X
!>          If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= n
!>          If INCX < 0, X(1-(n-i)*INCX) = x(i) for 1 <= i <= n
!>          If INCX = 0, x isn't a vector so there is no need to call
!>          this subroutine.  If you call it anyway, it will count x(1)
!>          in the vector norm N times.
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Edward Anderson, Lockheed Martin
!
!> \date August 2016
!
!> \ingroup nrm2
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
!>  Blue, James L. (1978)
!>  A Portable Fortran Program to Find the Euclidean Norm of a Vector
!>  ACM Trans Math Softw 4:15--23
!>  https://doi.org/10.1145/355769.355771
!>
!> \endverbatim
!>
!  =====================================================================
function SNRM2( n, x, incx ) 
   integer, parameter :: wp = kind(1.e0)
   real(wp) :: SNRM2
!
!  -- Reference BLAS level1 routine (version 3.9.1) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     March 2021
!
!  .. Constants ..
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: one  = 1.0_wp
   real(wp), parameter :: maxN = huge(0.0_wp)
!  ..
!  .. Blue's scaling constants ..
   real(wp), parameter :: tsml = real(radix(0._wp), wp)**ceiling( &
       (minexponent(0._wp) - 1) * 0.5_wp)
   real(wp), parameter :: tbig = real(radix(0._wp), wp)**floor( &
       (maxexponent(0._wp) - digits(0._wp) + 1) * 0.5_wp)
   real(wp), parameter :: ssml = real(radix(0._wp), wp)**( - floor( &
       (minexponent(0._wp) - digits(0._wp)) * 0.5_wp))
   real(wp), parameter :: sbig = real(radix(0._wp), wp)**( - ceiling( &
       (maxexponent(0._wp) + digits(0._wp) - 1) * 0.5_wp))
!  ..
!  .. Scalar Arguments ..
   integer :: incx, n
!  ..
!  .. Array Arguments ..
   real(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, ix
   logical :: notbig
   real(wp) :: abig, amed, asml, ax, scl, sumsq, ymax, ymin
!
!  Quick return if possible
!
   SNRM2 = zero
   if( n <= 0 ) return
!
   scl = one
   sumsq = zero
!
!  Compute the sum of squares in 3 accumulators:
!     abig -- sums of squares scaled down to avoid overflow
!     asml -- sums of squares scaled up to avoid underflow
!     amed -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     tbig -- values bigger than this are scaled down by sbig
!     tsml -- values smaller than this are scaled up by ssml
!
   notbig = .true.
   asml = zero
   amed = zero
   abig = zero
   ix = 1
   if( incx < 0 ) ix = 1 - (n-1)*incx
   do i = 1, n
      ax = abs(x(ix))
      if (ax > tbig) then
         abig = abig + (ax*sbig)**2
         notbig = .false.
      else if (ax < tsml) then
         if (notbig) asml = asml + (ax*ssml)**2
      else
         amed = amed + ax**2
      end if
      ix = ix + incx
   end do
!
!  Combine abig and amed or amed and asml if more than one
!  accumulator was used.
!
   if (abig > zero) then
!
!     Combine abig and amed if abig > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         abig = abig + (amed*sbig)*sbig
      end if
      scl = one / sbig
      sumsq = abig
   else if (asml > zero) then
!
!     Combine amed and asml if asml > 0.
!
      if ( (amed > zero) .or. (amed > maxN) .or. (amed /= amed) ) then
         amed = sqrt(amed)
         asml = sqrt(asml) / ssml
         if (asml > amed) then
            ymin = amed
            ymax = asml
         else
            ymin = asml
            ymax = amed
         end if
         scl = one
         sumsq = ymax**2*( one + (ymin/ymax)**2 )
      else
         scl = one / ssml
         sumsq = asml
      end if
   else
!
!     Otherwise all values are mid-range
!
      scl = one
      sumsq = amed
   end if
   SNRM2 = scl*sqrt( sumsq )
   return
end function
