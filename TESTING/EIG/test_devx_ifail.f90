! Check packed failure indices against the actual direct STEBZ/STEIN result.
! The split small block can fail inverse iteration with a coarse ABSTOL.
! Do not require nonconvergence on every implementation: any reported failures
! must still identify the corresponding columns after the drivers sort W/Z.
program test_devx_ifail
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: wp = kind(1.0d0), maxn = 7
  integer :: which, n, upper, mref, iref, family, calls, failures, witnesses
  real(wp) :: d(maxn), e(maxn), abstol, wref(maxn), sorted(maxn)
  logical :: expected(maxn)
  external :: dstebz, dstein, dstevx, dsyevx

  calls = 0
  failures = 0
  witnesses = 0
  do which = 1, 4
    call make_case(which, n, d, e, abstol, upper)
    call reference(n, d, e, abstol, upper, mref, iref, wref, sorted, expected)
    if (which == 3 .and. iref /= 0) error stop 2
    do family = 1, 3
      call driver(family, which, n, d, e, abstol, upper, mref, iref, sorted, expected)
    end do
  end do
  print '(a,i0,a,i0,a,i0)', 'EVX_INDEX_COMPLETE calls=', calls, &
       ' mismatches=', failures, ' genuine_failures=', witnesses
  if (failures /= 0) stop 1
contains
  subroutine make_case(which, n, d, e, abstol, upper)
    integer, intent(in) :: which
    integer, intent(out) :: n, upper
    real(wp), intent(out) :: d(maxn), e(maxn), abstol
    real(wp) :: factor
    integer :: j
    d = 0.0_wp
    e = 0.0_wp
    factor = 1.0_wp
    if (which /= 3) factor = scale(1.0_wp, -300)
    n = 5
    do j = 1, n
      d(j) = (1.3_wp + 0.37_wp * real(j-1, wp)) * factor
      if (j < n) e(j) = (0.27_wp / real(j, wp)) * factor
    end do
    abstol = 0.01_wp * factor
    if (which == 1) then
      d(3:7) = d(1:5)
      e(3:6) = e(1:4)
      d(1:2) = [10.0_wp, 12.0_wp]
      e(1:2) = 0.0_wp
      n = 7
    else if (which == 2) then
      d(6:7) = [-12.0_wp, -10.0_wp]
      e(5:6) = 0.0_wp
      n = 7
    end if
    upper = n
    if (n == 7) upper = n-1
  end subroutine

  subroutine reference(n, d, e, abstol, upper, m, info, w, ordered, expected)
    integer, intent(in) :: n, upper
    real(wp), intent(in) :: d(maxn), e(maxn), abstol
    integer, intent(out) :: m, info
    real(wp), intent(out) :: w(maxn), ordered(maxn)
    logical, intent(out) :: expected(maxn)
    real(wp) :: z(maxn,maxn), work(8*maxn)
    integer :: iblock(maxn), isplit(maxn), iwork(5*maxn), ifail(maxn)
    integer :: split, eigen_info, j, k, rank(maxn)
    w = 0.0_wp
    ordered = 0.0_wp
    z = 0.0_wp
    ifail = 0
    expected = .false.
    call dstebz('I', 'B', n, 0.0_wp, 0.0_wp, 1, upper, abstol, &
                d, e, m, split, w, iblock, isplit, work, iwork, eigen_info)
    if (eigen_info /= 0 .or. m /= upper) error stop 2
    call dstein(n, d, e, m, w, iblock, isplit, z, maxn, work, iwork, ifail, info)
    if (info < 0 .or. info > m) error stop 2
    ! Independent stable ranks; no copy of the drivers' selection-sort loop.
    do j = 1, m
      rank(j) = 1
      do k = 1, m
        if (w(k) < w(j) .or. (w(k) == w(j) .and. k < j)) rank(j) = rank(j)+1
      end do
      ordered(rank(j)) = w(j)
    end do
    do j = 1, info
      if (ifail(j) < 1 .or. ifail(j) > m) error stop 2
      if (expected(rank(ifail(j)))) error stop 2
      expected(rank(ifail(j))) = .true.
    end do
  end subroutine

  subroutine driver(family, which, n, d, e, abstol, upper, mref, iref, ordered, expected)
    integer, intent(in) :: family, which, n, upper, mref, iref
    real(wp), intent(in) :: d(maxn), e(maxn), abstol, ordered(maxn)
    logical, intent(in) :: expected(maxn)
    real(wp) :: dc(maxn), ec(maxn), a(maxn,maxn), w(maxn), z(maxn,maxn), work(8*maxn)
    integer :: iwork(5*maxn), ifail(maxn), info, m, j
    logical :: reported(maxn), valid
    character :: uplo
    character(8) :: name
    dc = d
    ec = e
    a = 0.0_wp
    z = 0.0_wp
    w = 0.0_wp
    ifail = 0
    if (family == 1) then
      name = 'DSTEVX'
      call dstevx('V', 'I', n, dc, ec, 0.0_wp, 0.0_wp, 1, upper, abstol, &
                  m, w, z, maxn, work, iwork, ifail, info)
    else
      do j = 1, n
        a(j,j) = d(j)
        if (j < n) then
          a(j,j+1) = e(j)
          a(j+1,j) = e(j)
        end if
      end do
      uplo = 'U'
      name = 'DSYEVXU'
      if (family == 3) then
        uplo = 'L'
        name = 'DSYEVXL'
      end if
      call dsyevx('V', 'I', uplo, n, a, maxn, 0.0_wp, 0.0_wp, 1, upper, abstol, &
                  m, w, z, maxn, work, size(work), iwork, ifail, info)
    end if
    calls = calls+1
    if (m /= mref .or. info /= iref) then
      failures = failures+1
      print '(a,a,i0)', trim(name), ' status differs from direct reference, case=', which
      return
    end if
    if (which <= 2 .and. info > 0) witnesses = witnesses+1
    if (.not. all(ieee_is_finite(w(1:m))) .or. .not. all(ieee_is_finite(z(1:n,1:m)))) then
      failures = failures+1
      return
    end if
    if (any(abs(w(1:m)-ordered(1:m)) > &
            8*epsilon(1.0_wp)*max(tiny(1.0_wp),abs(ordered(1:m))))) then
      failures = failures+1
      return
    end if
    reported = .false.
    valid = .true.
    do j = 1, info
      if (ifail(j) < 1 .or. ifail(j) > m) then
        valid = .false.
      else if (reported(ifail(j))) then
        valid = .false.
      else
        reported(ifail(j)) = .true.
      end if
    end do
    if (info == 0 .and. any(ifail(1:m) /= 0)) valid = .false.
    ! Only the first INFO entries are defined for a partial failure.
    if (.not. valid .or. any(reported .neqv. expected)) then
      failures = failures+1
      print '(a,a,i0,a,i0,a,7(i0,1x))', trim(name), ' IFAIL mismatch case=', which, &
           ' info=', info, ' indices=', ifail(1:m)
    end if
  end subroutine
end program
