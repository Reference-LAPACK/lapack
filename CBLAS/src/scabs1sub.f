c     scabs1.f
c
c     The program is a fortran wrapper for scabs1.
c
      subroutine scabs1sub(z, cabs1)
c
      external scabs1
      complex z
      real scabs1, cabs1
c
      cabs1=scabs1(z)
      return
      end
