subroutine perm_ice (f_ghz, t_k, perm_re, perm_im)

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

! Calculates the permittivity of ice
!*  Hufford (1991), Mishima et al. (1983) modified by Maetzler and Wegmueller (1987)

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

!*local variables
real (kind=jprb) :: tk, tc, theta, alpha, fac, betam, dbeta, beta

tk = min (t_k, 273.16_jprb)
tc = tk - 273.16_jprb

perm_re = 3.1884_jprb + 9.1e-04_jprb * tc

theta = 300.0_jprb / tk
alpha = 1.0e-04_jprb * (50.4_jprb + 62.0_jprb * (theta - 1.0_jprb)) * exp (-22.1_jprb * (theta - 1.0_jprb))

fac   = exp (335.0_jprb / tk)
betam = 0.0207_jprb * fac / (tk * (fac - 1.0_jprb) * (fac - 1.0_jprb)) + 1.16e-11_jprb * f_ghz  * f_ghz
dbeta = exp (-10.02_jprb + 0.0364_jprb * tc)
beta  = betam + dbeta

perm_im = alpha / f_ghz + beta * f_ghz

return
end subroutine perm_ice
