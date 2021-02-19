module LA_CONSTANTS
!
!  -- BLAS/LAPACK module --
!     May 06, 2016
!
!  Standard constants
!
   integer, parameter :: wp =    8
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
   character*1, parameter :: sprefix = 'D'
   character*1, parameter :: cprefix = 'Z'
!
!  Model parameters
!
   real(wp), parameter :: eps =  0.11102230246251565404E-015_wp
   real(wp), parameter :: ulp =  0.22204460492503130808E-015_wp
   real(wp), parameter :: safmin =  0.22250738585072013831E-307_wp
   real(wp), parameter :: safmax =  0.44942328371557897693E+308_wp
   real(wp), parameter :: smlnum =  0.10020841800044863890E-291_wp
   real(wp), parameter :: bignum =  0.99792015476735990583E+292_wp
   real(wp), parameter :: rtmin =  0.10010415475915504622E-145_wp
   real(wp), parameter :: rtmax =  0.99895953610111751404E+146_wp
!
!  Blue's scaling constants
!
   real(wp), parameter :: tsml =  0.14916681462400413487E-153_wp
   real(wp), parameter :: tbig =  0.19979190722022350281E+147_wp
   real(wp), parameter :: ssml =  0.44989137945431963828E+162_wp
   real(wp), parameter :: sbig =  0.11113793747425387417E-161_wp
end module LA_CONSTANTS
