module LA_CONSTANTS32
!
!  -- BLAS/LAPACK module --
!     May 06, 2016
!
!  Standard constants
!
   integer, parameter :: wp =    4
   real(wp), parameter :: zero = 0.0_wp
   real(wp), parameter :: half = 0.5_wp
   real(wp), parameter :: one = 1.0_wp
   real(wp), parameter :: two = 2.0_wp
   real(wp), parameter :: three = 3.0_wp
   real(wp), parameter :: four = 4.0_wp
   real(wp), parameter :: eight = 8.0_wp
   real(wp), parameter :: ten = 10.0_wp
   complex(wp), parameter :: czero = ( 0.0_wp, 0.0_wp )
   complex(wp), parameter :: chalf = ( 0.5_wp, 0.0_wp )
   complex(wp), parameter :: cone = ( 1.0_wp, 0.0_wp )
   character*1, parameter :: sprefix = 'S'
   character*1, parameter :: cprefix = 'C'
!
!  Model parameters
!
   real(wp), parameter :: eps =  0.5960464478E-07_wp
   real(wp), parameter :: ulp =  0.1192092896E-06_wp
   real(wp), parameter :: safmin =  0.1175494351E-37_wp
   real(wp), parameter :: safmax =  0.8507059173E+38_wp
   real(wp), parameter :: smlnum =  0.9860761315E-31_wp
   real(wp), parameter :: bignum =  0.1014120480E+32_wp
   real(wp), parameter :: rtmin =  0.3140184864E-15_wp
   real(wp), parameter :: rtmax =  0.3184525782E+16_wp
!
!  Blue's scaling constants
!
   real(wp), parameter :: tsml =  0.1084202172E-18_wp
   real(wp), parameter :: tbig =  0.4503599627E+16_wp
   real(wp), parameter :: ssml =  0.3777893186E+23_wp
   real(wp), parameter :: sbig =  0.1323488980E-22_wp
end module LA_CONSTANTS32
