function ice_density (i_type, f_ghz)

! This software was developed within the context of
! the EUMETSAT Satellite Application Facility on
! Numerical Weather Prediction (NWP SAF), under the
! Cooperation Agreement dated 25 November 1998, between
! EUMETSAT and the Met Office, UK, by one or more partners
! within the NWP SAF. The partners in the NWP SAF are
! the Met Office, ECMWF, KNMI and MeteoFrance.
!
! Copyright 2002, EUMETSAT, All Rights Reserved.

! This function implements the snow, graupel and ice density parametrization 
! from table II, Surussavadee and Staelin (2006), IEEE TGARS, 44, 2667-2678. 
! This gives a "density" as a function of frequency which is designed to 
! match scattering cross sections from Mie sphere computations to those from 
! the discrete-dipole approximation. This function can be called for all 
! hydrometeor types, and will just return a constant density for rain, clw etc.

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           10/02/2011  initial version (Alan Geer)

use parkind1, only: jprb, jpim
use mod_mie,  only: ice_a, ice_b

implicit none

! Interface
integer (kind=jpim), intent ( in) :: i_type
real    (kind=jprb), intent ( in) :: f_ghz

! Function definition
real (kind=jprb) :: ice_density ! [g cm^-3]

! Local variables
real (kind=jprb) :: f_thz

f_thz = f_ghz / 1000.0_JPRB
ice_density = ice_a(i_type) + f_thz * ice_b(i_type) 

return
end function ice_density
