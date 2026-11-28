! Regression test for the NaN pivot guard in SPBSTF and for the INFO = N+i
! that SSBGV, SSBGVD and SSBGVX document for a B that is not positive
! definite.
!
! A NaN is not .LE. 0, so SPBSTF used to factor a band matrix that holds a
! NaN into a NaN factor and return INFO = 0, and the generalized drivers
! went on to report a failure of the eigensolver (or, for SSBGVX, nothing at
! all) instead of INFO = N+i.  The xeigtstd drivers only generate positive
! definite B, so this program puts a NaN into every position of the band of
! a positive definite B and requires a nonzero INFO from SPBSTF (INFO = j for
! a NaN on the diagonal at j) and INFO = N+j from the drivers.
program test_spbstf_nan
  use, intrinsic :: iso_fortran_env, only: sp => real32
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  integer, parameter :: ncase = 6
  integer, parameter :: nn(ncase) = (/ 12, 11, 7, 6, 2, 1 /)
  integer, parameter :: kk(ncase) = (/ 3, 4, 1, 5, 1, 0 /)
  external :: spbstf, ssbgv, ssbgvd, ssbgvx
  character :: uplo
  integer :: icase, iuplo, n, kd, i, j, idiag, failures
  real(sp) :: nan

  nan = ieee_value(1.0_sp, ieee_quiet_nan)
  failures = 0

  do icase = 1, ncase
    n = nn(icase)
    kd = kk(icase)
    do iuplo = 1, 2
      uplo = merge('U', 'L', iuplo == 1)
      idiag = merge(kd + 1, 1, iuplo == 1)
      call finite_case()
      ! Every entry of the band, as (row of AB, column of AB).
      do j = 1, n
        do i = 1, kd + 1
          if (stored(i, j)) call nan_case(i, j, i == idiag)
        end do
      end do
    end do
  end do

  if (failures /= 0) then
    print '(a,i0,a)', 'test_spbstf_nan: ', failures, ' failures'
    stop 1
  end if
  print '(a)', 'test_spbstf_nan: all tests passed'

contains

  ! True if AB(i,j) is inside the band of the N-by-N matrix.
  logical function stored(i, j)
    integer, intent(in) :: i, j
    if (uplo == 'U') then
      stored = (i + j - kd - 1 >= 1)
    else
      stored = (i + j - 1 <= n)
    end if
  end function

  ! Positive definite B in band storage: 4 on the diagonal, 1 next to it.
  subroutine make_b(bb)
    real(sp), intent(out) :: bb(kd+1, n)
    bb = 0.0_sp
    if (uplo == 'U') then
      bb(kd+1, :) = 4.0_sp
      if (kd > 0) bb(kd, :) = 1.0_sp
    else
      bb(1, :) = 4.0_sp
      if (kd > 0) bb(2, :) = 1.0_sp
    end if
  end subroutine

  subroutine report(label, got, want)
    character(*), intent(in) :: label
    integer, intent(in) :: got, want
    print '(a,a,a,i0,a,i0,a,a,a,i0,a,i0)', 'FAIL ', label, ': N=', n, &
      ' KD=', kd, ' UPLO=', uplo, ' INFO=', got, ' expected ', want
    failures = failures + 1
  end subroutine

  ! Finite input: everything must succeed.
  subroutine finite_case()
    real(sp) :: ab(kd+1, n), bb(kd+1, n)
    integer :: info
    call make_b(bb)
    call spbstf(uplo, n, kd, bb, kd + 1, info)
    if (info /= 0) call report('SPBSTF finite', info, 0)
    call make_b(ab)
    call make_b(bb)
    call run_drivers(ab, bb, 0, 'finite')
  end subroutine

  ! NaN at AB(i,j) of B.  A NaN on the diagonal is the pivot itself, so
  ! SPBSTF must stop exactly there.  A NaN off the diagonal reaches the
  ! pivots of the rows and columns it is combined with through the rank-one
  ! updates, so SPBSTF must report some column.  The drivers add N to the
  ! column that SPBSTF reports; their INFO is compared only for a NaN on the
  ! diagonal, where that column is known.
  subroutine nan_case(i, j, diag)
    integer, intent(in) :: i, j
    logical, intent(in) :: diag
    real(sp) :: ab(kd+1, n), bb(kd+1, n)
    character(40) :: what
    integer :: info
    write (what, '(a,i0,a,i0,a)') 'NaN at AB(', i, ',', j, ')'
    call make_b(bb)
    bb(i, j) = nan
    call spbstf(uplo, n, kd, bb, kd + 1, info)
    if (diag) then
      if (info /= j) call report('SPBSTF ' // trim(what), info, j)
    else
      if (info < 1 .or. info > n) &
        call report('SPBSTF ' // trim(what), info, 1)
    end if
    if (diag) then
      call make_b(ab)
      call make_b(bb)
      bb(i, j) = nan
      call run_drivers(ab, bb, n + j, trim(what))
    end if
  end subroutine

  subroutine run_drivers(ab0, bb0, want, what)
    real(sp), intent(in) :: ab0(kd+1, n), bb0(kd+1, n)
    integer, intent(in) :: want
    character(*), intent(in) :: what
    real(sp) :: ab(kd+1, n), bb(kd+1, n), w(n), z(n, n), q(n, n)
    real(sp) :: work(7*n + 2*n*n + 1)
    integer :: iwork(5*n + 3), ifail(n), m, info, ivar
    character :: jobz
    do ivar = 1, 6
      ab = 0.5_sp * ab0
      bb = bb0
      jobz = merge('V', 'N', mod(ivar, 2) == 1)
      select case ((ivar + 1) / 2)
      case (1)
        call ssbgv(jobz, uplo, n, kd, kd, ab, kd + 1, bb, kd + 1, w, z, n, &
                   work, info)
        if (info /= want) call report('SSBGV ' // what // ' ' // jobz, info, &
                                 want)
      case (2)
        call ssbgvd(jobz, uplo, n, kd, kd, ab, kd + 1, bb, kd + 1, w, z, n, &
                    work, size(work), iwork, size(iwork), info)
        if (info /= want) call report('SSBGVD ' // what // ' ' // jobz, info, &
                                 want)
      case (3)
        call ssbgvx(jobz, 'A', uplo, n, kd, kd, ab, kd + 1, bb, kd + 1, q, &
                    n, 0.0_sp, 0.0_sp, 1, n, 0.0_sp, m, w, z, n, work, &
                    iwork, ifail, info)
        if (info /= want) call report('SSBGVX ' // what // ' ' // jobz, info, &
                                 want)
      end select
    end do
  end subroutine

end program test_spbstf_nan
