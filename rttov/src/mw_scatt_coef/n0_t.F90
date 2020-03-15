SUBroutine n0_t (i_type, lwc, temp, n0_mp)

! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------


use parkind1, only: jprb, jpim
!INTF_OFF
use mod_mie , only: t_fl, lwc2rr_a, lwc2rr_b, i_rain, i_snow, i_graupel
!INTF_ON
implicit none

!* common variables
integer (kind=jpim), intent ( in) :: i_type
real (kind=jprb), intent ( in) :: lwc, temp
real (kind=jprb), intent (out) :: n0_mp

!INTF_END

! local variables
real (kind=jprb)            :: rr

! End of header
!---------------------------------

rr = lwc2rr_a * lwc ** lwc2rr_b

! Parameterization of n0 vs T from Panegrossi et al. (1998), appendix.

if     (i_type == i_rain) then   

  n0_mp = 0.07_jprb * rr ** 0.37_jprb

elseif (i_type == i_snow) then  

  if (temp <= t_fl) then
    n0_mp = 0.04_jprb / exp (3.8E-02_jprb * (temp - t_fl))
  else
    n0_mp = 0.04_jprb / exp (8.8E-02_jprb * (temp - t_fl))
  endif

elseif (i_type == i_graupel) then   

  if (temp <= t_fl) then
    n0_mp = 0.04_jprb / exp (2.5E-02_jprb * (temp - t_fl))
  else
    n0_mp = 6.52E-06_jprb * (rr ** 0.872_jprb) / exp (7.5E-02_jprb * (temp - t_fl))
  endif

endif

return
end subroutine n0_t

