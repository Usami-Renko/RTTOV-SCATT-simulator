subroutine perm_melt (fv_melt, fm_melt, density_snow, &
                      perm_wat, perm_froz, alpha_melt, perm_core, perm_coat)

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


use parkind1, only: jprb
!INTF_OFF
use mod_mie, only : density_wat, density_air, density_ice, density_lmt, perm_air
!INTF_ON
!* Coated melting particle (Fabry and Szyrmer, JAS 1999)

implicit none

!* common variables
real    (kind=jprb), intent ( in) :: fv_melt, fm_melt, density_snow
real    (kind=jprb), intent (out) :: alpha_melt
complex (kind=jprb), intent ( in) :: perm_wat , perm_froz
complex (kind=jprb), intent (out) :: perm_core, perm_coat

!INTF_END

#include "vol_fracs.interface"
#include "mg_ellips.interface"

!* local variables
real (kind=jprb) :: alpha, beta, gamma, slope, offst
real (kind=jprb) :: density_snow_core, density_snow_coat, density_melt_core, density_melt_coat
real (kind=jprb) :: fv_ice_core, fv_air_core, fv_ice_coat, fv_air_coat

complex (kind=jprb) :: perm_ice_wat_core, perm_ice_wat_coat

!* eqs. 4, 5      
alpha =  0.5_jprb
beta  = -1.0_jprb
gamma = beta + 3.0_jprb     

density_snow_core = density_snow * (alpha ** beta)
density_snow_coat = density_snow * ((1.0_jprb - alpha ** gamma) / (1.0_jprb - alpha ** 3.0_jprb))

if (density_snow_core > density_ice) density_snow_core = density_ice
if (density_snow_coat > density_ice) density_snow_coat = density_ice
if (density_snow_core < density_lmt) density_snow_core = density_lmt
if (density_snow_coat < density_lmt) density_snow_coat = density_lmt

!* eqs. 6, 7      
density_melt_core = density_snow_core * density_wat / (fm_melt * density_snow_core + (1.0_jprb - fm_melt) * density_wat)
density_melt_coat = density_snow_coat * density_wat / (fm_melt * density_snow_coat + (1.0_jprb - fm_melt) * density_wat) 

slope = 0.5_jprb / (1.0_jprb - density_snow_coat / density_snow_core)
offst = 1.0_jprb - slope

!* eq. 8
alpha_melt = slope * density_melt_coat / density_melt_core + offst 

!* reduction of ice/air volume fractions with melting
call vol_fracs (density_ice, density_air, density_snow_core, fv_melt, fv_ice_core, fv_air_core)
call vol_fracs (density_ice, density_air, density_snow_coat, fv_melt, fv_ice_coat, fv_air_coat)     
              
!* inclusions: relative fv's
call mg_ellips (fv_ice_core / (fv_ice_core + fv_melt), perm_wat, perm_froz, perm_ice_wat_core)
call mg_ellips (fv_ice_coat / (fv_ice_coat + fv_melt), perm_wat, perm_froz, perm_ice_wat_coat)
 
!* core: air in (ice in water); coat: (ice in water) in air         
call mg_ellips (           fv_air_core, perm_ice_wat_core, perm_air         , perm_core)         
call mg_ellips (1.0_jprb - fv_air_core, perm_air         , perm_ice_wat_coat, perm_coat)

return
end subroutine perm_melt
