! SPDX-FileCopyrightText: The LAPACK Authors
! SPDX-License-Identifier: BSD-3-Clause
!
! Regression tests for the RQ solve and residual test helpers.
program ztest_rq_helpers
  use, intrinsic :: ieee_arithmetic
  implicit none
  integer, parameter :: wp = kind(1.0d0)
  complex(wp) :: a(2,3), b(3,2), x(3,2), tau(2), work(128)
  real(wp) :: large, small, resid, rwork(3), bad(2)
  integer :: info, i, j, k
  external :: zgerqs, zget02
  character :: trans(3) = ['N', 'T', 'C']

  ! R has both ends of the normal range in one right-hand side.
  ! Uniform scaling of B would lose the first solution component.
  large = 0.75_wp*huge(1.0_wp)
  small = tiny(1.0_wp)
  a = 0
  a(1,2) = small
  a(2,3) = large
  tau = 0
  b = 0
  b(2,1) = small
  b(3,1) = cmplx(-0.95_wp*large, 0.95_wp*large, wp)
  b(2,2) = 2*small
  b(3,2) = 0.5_wp*large
  x = 0
  x(2,1) = 1
  x(3,1) = cmplx(-0.95_wp, 0.95_wp, wp)
  x(2,2) = 2
  x(3,2) = 0.5_wp
  call zgerqs(2, 3, 2, a, 2, tau, b, 3, work, 128, info)
  if (info /= 0) stop 1
  if (.not. all(abs(b-x) <= 8*epsilon(1.0_wp))) stop 2

  ! Check ordinary magnitudes as well.
  a(1,2) = 2
  a(2,3) = 4
  b = 0
  b(2,:) = 2*x(2,:)
  b(3,:) = 4*x(3,:)
  call zgerqs(2, 3, 2, a, 2, tau, b, 3, work, 128, info)
  if (info /= 0) stop 3
  if (.not. all(abs(b-x) <= 8*epsilon(1.0_wp))) stop 4

  bad = [ieee_value(0.0_wp, ieee_quiet_nan), &
         ieee_value(0.0_wp, ieee_positive_inf)]
  do k = 1, 3
    do j = 1, 3
      do i = 1, 2
        a = 1
        x = 1
        b = 1
        if (j == 1) a(1,1) = bad(i)
        if (j == 2) x(1,1) = bad(i)
        if (j == 3) b(1,1) = bad(i)
        call zget02(trans(k), 1, 1, 1, a, 2, x, 3, b, 3, rwork, resid)
        if (.not. (resid > 30)) stop 5
      end do
    end do
    a = 1
    x = 1
    b = 1
    call zget02(trans(k), 1, 1, 1, a, 2, x, 3, b, 3, rwork, resid)
    if (resid /= 0) stop 6
  end do
end program
