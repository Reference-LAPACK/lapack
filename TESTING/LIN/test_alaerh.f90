! SPDX-FileCopyrightText: The LAPACK Authors
! SPDX-License-Identifier: BSD-3-Clause
!
! Unexpected success must be diagnosed like any other wrong INFO.
program test_alaerh
  implicit none
  character :: precision(4) = ['S', 'D', 'C', 'Z']
  character(2) :: storage(3) = ['PB', 'PP', 'PT']
  character(3) :: path
  character(6) :: subnam
  integer :: p, s

  do p = 1, size(precision)
    do s = 1, size(storage)
      path = precision(p)//storage(s)
      subnam = path//'TRF'
      call check(0, 3, 0)
      call check(2, 3, 1)
      call check(2, 0, 1)
      call check(2, 2, 1)
      call check(0, 0, 0)
    end do
  end do

contains

  subroutine check(info, infoe, prior_errors)
    integer, intent(in) :: info, infoe, prior_errors
    integer :: unit, nerrs, nfail, status, lines, diagnostics
    character(256) :: line
    character(5) :: actual
    character(2) :: expected
    external :: alaerh

    nfail = 0
    nerrs = prior_errors
    open(newunit=unit, status='scratch', action='readwrite')
    call alaerh(path, subnam, info, infoe, 'U', 5, 5, 1, 1, 1, &
                9, nfail, nerrs, unit)
    if (nfail /= 0) stop 1
    rewind(unit)
    lines = 0
    diagnostics = 0
    write(actual, '(I5)') info
    write(expected, '(I2)') infoe
    do
      read(unit, '(A)', iostat=status) line
      if (status < 0) exit
      if (status /= 0) stop 2
      lines = lines + 1
      if (index(line, subnam) == 0) cycle
      if (index(line, ' *** ') == 0) cycle
      diagnostics = diagnostics + 1
      if (index(line, '='//actual) == 0) stop 3
      if (info /= infoe .and. infoe /= 0) then
        if (index(line, 'instead of '//expected) == 0) stop 4
      end if
    end do
    close(unit)
    if (info == 0 .and. infoe == 0) then
      if (nerrs /= prior_errors .or. lines /= 0) stop 5
    else
      if (nerrs /= prior_errors+1 .or. diagnostics /= 1) stop 6
    end if
  end subroutine
end program
