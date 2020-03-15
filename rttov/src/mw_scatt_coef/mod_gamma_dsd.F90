subroutine mod_gamma_dsd (lwc, diameter, nr)

! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2002, EUMETSAT, All Rights Reserved.

! calculates the modified gamma size distribution

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)

use parkind1, only: jprb

implicit none

!* common variables
real (kind=jprb), intent ( in) :: lwc, diameter
real (kind=jprb), intent (out) :: nr

!INTF_END

!* local variables	
real (kind=jprb) :: r_micron

nr       = 0.0_jprb
r_micron = diameter * 0.5e+04_jprb

nr = 7.676_jprb * (r_micron ** 2.0_jprb) / exp (0.425_jprb * r_micron)

nr = nr * 1.0e+04_jprb * lwc / 1.31_jprb

return
end subroutine mod_gamma_dsd
