subroutine gamma_dsd (i_type, temp, lwc, mue, density, n0_fac, ll_dsd, n0, d0, lambda, &
  & a, b, gamma1plusb)

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

! This program calculates the size distribution given the water content and density

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!           29/10/2009   hunting down differences with Alan's code, making changes to get it same. AMD
!           12/03/2010   restructuring (Alan Geer)
!           15/03/2013   fully-flexible PSD, density and shape options (Alan Geer)

use parkind1, only: jprb, jpim, jplm
!INTF_OFF
use mod_mie  , only: i_aggregate, not_available
!INTF_ON

implicit none

!* common variables
integer (kind=jpim), intent ( in) :: i_type
real    (kind=jprb), intent ( in) :: temp, lwc, mue, density(:), n0_fac
logical (kind=jplm), intent ( in) :: ll_dsd
real    (kind=jprb), intent (out) :: n0, d0, lambda
real    (kind=jprb), intent ( in) :: a, b, gamma1plusb

!INTF_END
#include "n0_t.interface"

!* local variables	
real        (kind=jprb) :: gamma_fct, fac, n0_mp
integer     (kind=jpim) :: imue, integer_limit

integer_limit = nint (mue) + 4 - 1

gamma_fct = 1.0_jprb
do imue = 1, integer_limit
  gamma_fct = gamma_fct * imue
end do

if (i_type /= i_aggregate .and. ll_dsd) then
  call n0_t (i_type, lwc, temp, n0_mp)
else
  n0_mp = 0.08_jprb
endif

n0 = n0_fac * n0_mp * exp (3.2_jprb * mue)

if(mue == 0._jprb) then

  ! Marshall-Palmer

  ! NB this maths is not appropriate for modified gamma distributions, just basic e.g. Marshall-Palmer
  ! (a,b are in SI.). 
  if(gamma1plusb == not_available .or. a == not_available) then
    write(*,*) 'Your density option is currently incompatible with the Marshall-Palmer PSD'
    stop
  endif

  ! Computation is in SI and then converted back to cm
  lambda = (((n0*1E8_JPRB)*a*gamma1plusb/(lwc/1E3_JPRB))**(1.0_JPRB/(1.0_JPRB+b)))/100.0_JPRB

else

  ! Full Gamma PSD

  if(any(density /= density(1))) then
    write(*,*) 'Your density option is currently incompatible with the Gamma PSD'
    stop
  endif

  !* pi * rho_x / (lwc)
  fac = 3.141592654_jprb * density(1) / (lwc * 6.0e-06_jprb)

  d0 = (((3.67_jprb + mue) ** (4.0_jprb + mue)) / (fac * gamma_fct * n0)) ** (1.0_jprb / (4.0_jprb + mue))
  lambda = (3.67_jprb + mue) / d0

endif

return
end subroutine gamma_dsd
