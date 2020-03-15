subroutine density_all (i_type, i_scat, particle_shape, dens, f_ghz, Dinmeters, &
  & density, a, b, gamma1plusb)

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


! Centralises all density-related computations. Outputs are (1) density
! itself and (2) the a/x and b/y coefficients of the mass-density relation.
! and other parameters for when the density is manipulated analytically. 
!
! A sharp observer will note that in some cases in this code the
! mass-density relation is specified in two separate ways, e.g. as
! an implicit definition of the x and y in m = xD^y in SI units, or as
! an explicit density definition such as density(:) = 0.132/D(:)/1000.0
! This disticntion is maintained for backward compatability, i.e. to 
! enable the exact re-generation of Mie tables created with earlier versions
! of the code (where density computations were scattered all about). 
! Ideally in the future, all density and mass properties will be calculated 
! directly from x and y if available.

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           03/2013     First version (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)

use parkind1, only: jprb, jpim
!INTF_OFF
use parkind1, only: jplm
use mod_mie, only:  density_liq, pi, dens_wilson_ballard_1999, dens_jones_1995, &
 & dens_brown_francis_1995, dens_fixed, dens_ss_2006, i_clw, i_rain, not_available, i_liu, i_arts
use mod_arts, only: alpha_arts, beta_arts
!INTF_ON

implicit none

integer (kind=jpim), intent ( in) :: i_type         ! Hydrometeor type
integer (kind=jpim), intent ( in) :: i_scat         ! Scattering computation type
integer (kind=jpim), intent ( in) :: particle_shape ! particle shape (e.g. Liu DDA habit number)
integer (kind=jpim), intent ( in) :: dens           ! Choice of density parametrisation
real    (kind=jprb), intent ( in) :: f_ghz          ! Frequency, [GHz]
real    (kind=jprb), intent ( in), optional :: Dinmeters(:)      ! Diameter or maximum dimension [m]
real    (kind=jprb), intent (out), optional :: density(:)        ! Density,   [g/m3]
real    (kind=jprb), intent (out), optional :: a, b, gamma1plusb ! mass = aD^b [SI]

!INTF_END

integer (kind=jpim) :: i_dia
real    (kind=jprb) :: x, y, gamma_local
logical (kind=jplm) :: lDD, lD

real    (kind=jprb) :: ice_density ! a function
#include "liu_density.interface"

lDD = present(Dinmeters) .and. present(density)
lD  = present(density) 

x = not_available
y = not_available
gamma_local = not_available

! The formulae are typically in SI units, translated to g/cm3 by dividing by 1000
if ( (i_scat == i_liu) .or. (i_scat == i_arts) ) then

  if (i_scat == i_liu) then
    ! Liu shapes (a,b are in SI.). See p. 27-139
    call liu_density(particle_shape, x, y, gamma_local)
  else
    x = alpha_arts(particle_shape)
    y = beta_arts(particle_shape)
  endif

  if( lDD ) then
    do i_dia = 1, size(Dinmeters)

      ! Coefficients are in SI, so convert dia_froz to metres
      density(i_dia) = 6.0_JPRB/pi*x*(Dinmeters(i_dia))**(y-3.0_JPRB)

    enddo 

    ! Then convert the density from kg m-3 to g cm-3 as is used more generally here
    density(:) = density(:) /1e3_JPRB
  endif

elseif (dens == dens_wilson_ballard_1999) then

  ! 1:density = 0.132*D-1, 
  if( lDD ) then
    density(:) = (0.132_jprb/Dinmeters(:))/1000.0_jprb
  endif

  ! set x and y to be values used in UM
  x=0.069
  y=2.0
  gamma_local = 3.0_jprb

elseif (dens == dens_jones_1995) then

  ! 2:density = 8.74E-4*exp(-0.625D2) + 4.5E-5
  if( lDD ) then
    density(:) = (0.000874_jprb*exp(-0.625_jprb*(Dinmeters(:)*Dinmeters(:))) + 0.000045_jprb)/1000.0_jprb
  endif

  ! x, y and gamma are not valid here!

elseif (dens == dens_brown_francis_1995) then

  ! 3:density = 0.035*D-1.1
  if( lDD ) then
    density(:) = (0.035_jprb/(Dinmeters(:)**1.1_jprb))/1000.0_jprb
  endif

  ! set x and y to be the values from Brown and Francis
  x = 0.0185_jprb
  y = 1.9
  gamma_local = 1.82736_jprb

elseif (dens == dens_fixed) then

  ! Fixed density (g/m3)
  if(lD) density(:) = density_liq(i_type)

  ! Fixed density sphere (x and y in SI units)
  x = 1e3_jprb * density_liq (i_type) * pi / 6.0_jprb
  y = 3.0_jprb
  gamma_local = 6.0_jprb

elseif (dens == dens_ss_2006) then

  ! Surrusavadee and Staelin (2006) ice factor (g/m3)
  if(lD) density(:) = ice_density(i_type, f_ghz)

  ! Fixed density sphere (x and y in SI units)
  x = 1e3_jprb * ice_density(i_type, f_ghz) * pi / 6.0_jprb
  y = 3.0
  gamma_local = 6.0_jprb

endif

if (present(a)) a=x
if (present(b)) b=y
if (present(gamma1plusb)) gamma1plusb = gamma_local 

! put in upper limit so density cannot be > 1.0
if(ld .and. i_type /= i_clw .and. i_type /= i_rain) then
  where (density > 0.92_jprb) density = 0.92_jprb
endif

return
end subroutine density_all
