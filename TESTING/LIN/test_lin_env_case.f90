! The test ILAENV must recognize short, padded and mixed-case QR/LQ names.
program test_lin_env_case
  implicit none
  integer, external :: ilaenv
  external :: xlaenv
  integer :: p,s,mode,n3,want,got,cases
  character, parameter :: prec(4)=['S','D','C','Z']
  character(6) :: name
  call xlaenv(1,17)
  call xlaenv(2,23)
  call xlaenv(3,31)
  call xlaenv(4,37)
  call xlaenv(5,41)
  cases=0
  do p=1,4
    do mode=1,6
      select case(mode)
      case(1)
        name=prec(p)//'GEQR '
      case(2)
        name=prec(p)//'geqr '
      case(3)
        name=prec(p)//'GeQr '
      case(4)
        name=prec(p)//'GELQ '
      case(5)
        name=prec(p)//'gelq '
      case(6)
        name=prec(p)//'GeLq '
      end select
      do n3=0,2
        want=17
        if(n3==2) want=23
        do s=1,5
          got=ilaenv(s,name,'N',10,5,n3,-1)
          if(got/=want) stop 1
          ! TRIM produces the five-character actual used by xGEQR/xGELQ.
          got=ilaenv(s,trim(name),'N',10,5,n3,-1)
          if(got/=want) stop 2
          cases=cases+2
        end do
      end do
    end do
  end do
  ! Other names retain the ordinary ISPEC-indexed dispatch.
  if(ilaenv(3,'DGEQRF','N',10,5,0,-1)/=31) stop 3
  ! Short names must not read beyond the actual argument.
  if(ilaenv(2,'D','N',10,5,0,-1)/=23) stop 4
  if(ilaenv(2,'','N',10,5,0,-1)/=23) stop 5
  print *, 'LIN environment name case: PASS',cases+3
end program
