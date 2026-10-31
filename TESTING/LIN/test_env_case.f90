! Environment queries must not depend on the case of routine names.
program test_env_case
  implicit none
  character, parameter :: prec(4) = ['S','D','C','Z']
  character(14), parameter :: suffixes(17) = [character(14) :: &
       'GEQRF','GEQR','GELQ','GEQP3RK','GEHRD','GEBRD','GETRF','TZRZF', &
       'ORMRZ','UNMRZ','SYTRD','HETRD','SYTRD_2STAGE','HETRD_2STAGE', &
       'LAORHR_GETRFNP','TREVC','TRSYL']
  character(6), parameter :: qrnames(5) = [character(6) :: &
       'HSEQR','LAQR0','GGHRD','GGHD3','TGEXC']
  character(12), parameter :: stages(6) = [character(12) :: &
       'SYTRD_2STAGE','HETRD_2STAGE','SYTRD_SY2SB','HETRD_HE2HB', &
       'GEBRD_GE2GB','GEBRD_GB2BD']
  integer, parameter :: sizes(9) = [1,14,16,64,150,256,590,3000,6000]
  character(16) :: name, variant
  integer :: p,j,k,s,mode,want,got,cases
  integer, external :: ilaenv, iparmq, iparam2stage
  cases = 0
  do p=1,4
    do j=1,size(suffixes)
      name=prec(p)//trim(suffixes(j))
      do k=1,size(sizes)
        do s=1,3
          want=ilaenv(s,trim(name),'N',sizes(k),sizes(k)/2+1,0,-1)
          do mode=1,2
            variant=case_name(name,mode)
            got=ilaenv(s,trim(variant),'N',sizes(k),sizes(k)/2+1,0,-1)
            call check(want,got,'ILAENV',name,variant,s)
          end do
        end do
      end do
    end do
    do j=1,size(qrnames)
      name=prec(p)//trim(qrnames(j))
      do k=1,size(sizes)
        do s=12,17
          want=iparmq(s,trim(name),'N',sizes(k),1,sizes(k),-1)
          do mode=1,2
            variant=case_name(name,mode)
            got=iparmq(s,trim(variant),'N',sizes(k),1,sizes(k),-1)
            call check(want,got,'IPARMQ',name,variant,s)
          end do
        end do
      end do
    end do
    do j=1,size(stages)
      name=prec(p)//trim(stages(j))
      do k=1,size(sizes)
        do s=17,21
          want=iparam2stage(s,trim(name),'N',sizes(k),32,16,0)
          do mode=1,2
            variant=case_name(name,mode)
            got=iparam2stage(s,trim(variant),'N',sizes(k),32,16,0)
            call check(want,got,'IPARAM2STAGE',name,variant,s)
          end do
        end do
      end do
    end do
  end do
  print *, 'Environment name case: PASS',cases
contains
  function case_name(input,mode) result(output)
    character(*),intent(in) :: input
    integer,intent(in) :: mode
    character(len(input)) :: output
    character(*),parameter :: upper='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(*),parameter :: lower='abcdefghijklmnopqrstuvwxyz'
    integer :: i,k
    output=input
    do i=1,len(input)
      k=index(upper,input(i:i))
      if(k>0 .and. (mode==1.or.mod(i,2)==0)) output(i:i)=lower(k:k)
    end do
  end function
  subroutine check(want,got,routine,name,variant,s)
    integer,intent(in) :: want,got,s
    character(*),intent(in) :: routine,name,variant
    if(want/=got) then
      print *, routine,': ',trim(name),' / ',trim(variant),s,want,got
      stop 1
    end if
    cases=cases+1
  end subroutine
end program
