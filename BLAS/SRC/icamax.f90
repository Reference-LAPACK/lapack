!> \brief \b ICAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ICAMAX(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>  ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
!>          X is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of X
!> \endverbatim
!
!  Authors:
!  ========
!
!> James Demmel, University of California Berkeley, USA
!> Weslley Pereira, National Renewable Energy Laboratory, USA
!
!> \ingroup iamax
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!> James Demmel et al. Proposed Consistent Exception Handling for the BLAS and
!> LAPACK, 2022 (https://arxiv.org/abs/2207.09281).
!>
!> \endverbatim
!>
!  =====================================================================
integer function icamax(n, x, incx)
   implicit none
   integer, parameter :: wp = kind(1.e0)
!
!  -- Reference BLAS level1 routine --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  .. Constants ..
   real(wp), parameter :: hugeval = huge(0.0_wp)
!
!  .. Scalar Arguments ..
   integer :: n, incx
!
!  .. Array Arguments ..
   complex(wp) :: x(*)
!  ..
!  .. Local Scalars ..
   integer :: i, j, ix, jx
   real(wp) :: val, smax
   logical :: scaledsmax
!  ..
!  .. Intrinsic Functions ..
   intrinsic :: abs, aimag, huge, real
!
!  Quick return if possible
!
   icamax = 0
   if (n < 1 .or. incx < 1) return
!
   icamax = 1
   if (n == 1) return
!
   icamax = 0
   scaledsmax = .false.
   smax = -1
!
!  scaledsmax = .true. indicates that x(icamax) is finite but
!  abs(real(x(icamax))) + abs(aimag(x(icamax))) overflows
!
   if (incx == 1) then
      ! code for increment equal to 1
      do i = 1, n
         if (x(i) /= x(i)) then
            ! return when first NaN found
            icamax = i
            return
         elseif (abs(real(x(i))) > hugeval .or. abs(aimag(x(i))) > hugeval) then
            ! keep looking for first NaN
            do j = i+1, n
               if (x(j) /= x(j)) then
                  ! return when first NaN found
                  icamax = j
                  return
               endif
            enddo
            ! record location of first Inf
            icamax = i
            return
         else ! still no Inf found yet
            if (.not. scaledsmax) then
               ! no abs(real(x(i))) + abs(aimag(x(i))) = Inf yet
               val = abs(real(x(i))) + abs(aimag(x(i)))
               if (val > hugeval) then
                  scaledsmax = .true.
                  smax = 0.25*abs(real(x(i))) + 0.25*abs(aimag(x(i)))
                  icamax = i
               elseif (val > smax) then ! everything finite so far
                  smax = val
                  icamax = i
               endif
            else ! scaledsmax
               val = 0.25*abs(real(x(i))) + 0.25*abs(aimag(x(i)))
               if (val > smax) then
                  smax = val
                  icamax = i
               endif
            endif
         endif
      end do
   else
      ! code for increment not equal to 1
      ix = 1
      do i = 1, n
         if (x(ix) /= x(ix)) then
            ! return when first NaN found
            icamax = i
            return
         elseif (abs(real(x(ix))) > hugeval .or. abs(aimag(x(ix))) > hugeval) then
            ! keep looking for first NaN
            jx = ix + incx
            do j = i+1, n
               if (x(jx) /= x(jx)) then
                  ! return when first NaN found
                  icamax = j
                  return
               endif
               jx = jx + incx
            enddo
            ! record location of first Inf
            icamax = i
            return
         else ! still no Inf found yet
            if (.not. scaledsmax) then
               ! no abs(real(x(ix))) + abs(aimag(x(ix))) = Inf yet
               val = abs(real(x(ix))) + abs(aimag(x(ix)))
               if (val > hugeval) then
                  scaledsmax = .true.
                  smax = 0.25*abs(real(x(ix))) + 0.25*abs(aimag(x(ix)))
                  icamax = i
               elseif (val > smax) then ! everything finite so far
                  smax = val
                  icamax = i
               endif
            else ! scaledsmax
               val = 0.25*abs(real(x(ix))) + 0.25*abs(aimag(x(ix)))
               if (val > smax) then
                  smax = val
                  icamax = i
               endif
            endif
         endif
         ix = ix + incx
      end do
   endif
end
