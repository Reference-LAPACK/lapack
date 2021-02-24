module LA_CONSTANTS32
!
!  -- BLAS/LAPACK module --
!     May 06, 2016
!
!  Standard constants
!
   real, parameter :: zero = 0.0
   real, parameter :: half = 0.5
   real, parameter :: one = 1.0
   real, parameter :: two = 2.0
   real, parameter :: three = 3.0
   real, parameter :: four = 4.0
   real, parameter :: eight = 8.0
   real, parameter :: ten = 10.0
   complex, parameter :: czero = ( 0.0, 0.0 )
   complex, parameter :: chalf = ( 0.5, 0.0 )
   complex, parameter :: cone = ( 1.0, 0.0 )
   character*1, parameter :: sprefix = 'S'
   character*1, parameter :: cprefix = 'C'
!
!  Model parameters
!
   real, parameter :: eps =  0.5960464478E-07
   real, parameter :: ulp =  0.1192092896E-06
   real, parameter :: safmin =  0.1175494351E-37
   real, parameter :: safmax =  0.8507059173E+38
   real, parameter :: smlnum =  0.9860761315E-31
   real, parameter :: bignum =  0.1014120480E+32
   real, parameter :: rtmin =  0.3140184864E-15
   real, parameter :: rtmax =  0.3184525782E+16
!
!  Blue's scaling constants
!
   real, parameter :: tsml =  0.1084202172E-18
   real, parameter :: tbig =  0.4503599627E+16
   real, parameter :: ssml =  0.3777893186E+23
   real, parameter :: sbig =  0.1323488980E-22
end module LA_CONSTANTS32
