subroutine perm_water (f_ghz, t_k, perm_re, perm_im)

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

! Calculates the permittivity of water
!* Liebe et al. (1989)

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)

use parkind1, only: jprb

implicit none

!* common variables
real (kind=jprb), intent ( in) :: f_ghz, t_k
real (kind=jprb), intent (out) :: perm_re, perm_im

!INTF_END

!* local variables
real (kind=jprb)  :: eps0, eps1, eps2, fp, fs, fac1, fac2, fac3, fac4, fac5, fac6
real (kind=jprb)  :: theta

theta = 300.0_jprb / t_k
fac1  = theta - 1.0_jprb
fac2  = fac1 * fac1

eps0 = 77.66_jprb + 103.3_jprb * fac1
eps1 = 5.48_jprb
eps2 = 3.51_jprb

fp = 20.09_jprb -  142.0_jprb * fac1 + 294.0_jprb * fac2
fs = 590.0_jprb - 1500.0_jprb * fac1

fac3 = f_ghz / fp
fac4 = 1.0_jprb + fac3 * fac3
fac5 = f_ghz / fs
fac6 = 1.0_jprb + fac5 * fac5

perm_re = (eps0 - eps1) / fac4 + (eps1 - eps2) / fac6 + eps2
perm_im = (eps0 - eps1) * fac3 / fac4 + (eps1 - eps2) * fac5 / fac6

return
end subroutine perm_water
