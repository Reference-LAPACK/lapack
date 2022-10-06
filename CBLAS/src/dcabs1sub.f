c     dcabs1.f
c
c     The program is a fortran wrapper for dcabs1.
c
      subroutine dcabs1sub(z, cabs1)
c
      external dcabs1
      double complex z
      double precision dcabs1, cabs1
c
      cabs1=dcabs1(z)
      return
      end
