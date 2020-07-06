function permittivity (i_type, f_ghz, temp, ipermwat, density)

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

! Current Code Owner: SAF NWP

! Inputs
!   i_type: hydrometeor type as in mod_mie.F90
!   f_ghz:  frequency [GHz]
!   temp:   temperature [K]
!   ipermwat: (OPTIONAL) water permittivity option as in mod_mie.F90
!   density:  (OPTIONAL) if present, and i_type is a frozen hydrometeor, 
!                        return effective permittivity of ice-air mixture

! History:
! Version   Date        Comment
! -------   ----        -------
!           10/03/2010  extracted from create_tables_spectra (Alan Geer)
!           03/03/2011  "ice-factor" density option
!           13/03/2013  Fully-flexible PSD, density and shape (Alan Geer)
!           10/01/2018  add option to select different liquid water permittivity models (Katrin Lonitz)
!           16/01/2018  Every permittivity calculation now goes through this routine, 
!                       inc. radar reflectivity factor and melting layer computations (Alan Geer)

use parkind1, only: jprb, jpim
!INTF_OFF
use mod_mie,  only: i_rain, i_clw, perm_air, density_air, density_ice, permwat_liebe_89, &
  & permwat_rosenkranz_15 , permwat_TKC_16
!INTF_ON

implicit none

! Interface
integer (kind=jpim), intent ( in) :: i_type
real    (kind=jprb), intent ( in) :: f_ghz
real    (kind=jprb), intent ( in) :: temp
integer (kind=jpim), optional, intent ( in) :: ipermwat
real    (kind=jprb), optional, intent ( in) :: density

! Function definition
complex (kind=jprb) :: permittivity

!INTF_END

! Local variables
real    (kind=jprb) :: perm_re, perm_im
real (kind=jprb)    :: fv_air, fv_ice
complex (kind=jprb) :: perm_froz

! Interfaces of called functions
#include "perm_water_liebe_89.interface"
#include "perm_water_rosenkranz_15.interface"
#include "perm_water_TKC_16.interface"
#include "perm_ice.interface"
#include "mg_ellips.interface"

!* permittivities
if (i_type == i_rain .or. i_type == i_clw) then

  if( .not. present(ipermwat) ) then
    stop 'Permittivity model of liquid water must be specified'
  endif

  !* Liquid
  select case(ipermwat)

    case(permwat_liebe_89)
      call perm_water_liebe_89 (f_ghz, temp, perm_re, perm_im)

    case(permwat_rosenkranz_15)
      call perm_water_rosenkranz_15 (f_ghz, temp, perm_re, perm_im)

    case(permwat_TKC_16)
      call perm_water_TKC_16 (f_ghz, temp, perm_re, perm_im)
      
    case default
      stop 'Unknown permittivity model of liquid water'

  end select
  permittivity = cmplx (perm_re, perm_im, jprb)

else

  !* Frozen
  call perm_ice (f_ghz, temp, perm_re, perm_im)
  perm_froz = cmplx (perm_re, perm_im, jprb)

  if( present(density) ) then
  
    ! Use mixing rule to compute effective permittivity of ice-air mixture
    fv_ice = (density - density_air) / (density_ice - density_air)
    fv_air = 1.0_jprb - fv_ice
    call mg_ellips (fv_air, perm_froz, perm_air, permittivity)
  
  else
  
    ! Return solid ice permittivity
    permittivity = perm_froz
    
  endif

endif

return
end function permittivity
