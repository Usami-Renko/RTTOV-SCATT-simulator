subroutine set_spectra (i_type, lwc, temp, f_ghz, dia_froz, nd, psd, ll_dsd, dens, regime, &
  & i_scat, particle_habit, ll_variable_d)

! Description:
!
!   Computes the particle size distribution appropriate to a given water content. Many 
!   options are available:
!
!      Field et al. (2005) - only for totalice (option 3)
!      Field et al. (2007) - only for totalice (option 4) or Liu DDA "snow"
!      Marshall-Palmer     - snow, rain, graupel, aggregate, totalice (option 2)
!      Modified Gamma      - CLW, CIW, totalice (option 1)
!
!   In all cases, the particle density relation (which may be constant or a function of 
!   particle size) is a basic input, and is configurable for totalice.
!
!   IN: i_type     - hydrometeor type (see mod_mie.F90)   
!       lwc        - liquid water content [g m^-3]
!       temp       - temperature          [K]
!       f_ghz      - frequency            [GHz]
!       dia_froz   - particle diameter or maximum dimension [cm]
!       dens       - density parametrization (totalice only)
!       psd        - PSD parametrization
!       regime     - Field et al. (2007) PSD regime
!       ll_dsd     - Panegrossi et al. n0 vs T
!       i_scat     - scattering computation type (see mod_mie.F90)
!       particle_habit - shape ID for Liu (2008) DDA shapes or Baran ensemble etc.
!         
!   OUT: nd        - size distribution [cm^-4]
!
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
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!           12/03/2010  Put all size distributions here (Alan Geer)
!           09/02/2011  Simplify normalisation. Documentation. (Alan Geer)
!           15/03/2013  Fully-flexible PSD, density and shape options (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)

use parkind1, only: jprb, jpim, jplm
use mod_mie,  only: n_dia
!INTF_OFF
! use mod_mie,  only: mue, density_liq, n0_fac, &
!                & onethird, gm32gcm3, pi, i_rain, i_snow, i_graupel, &
!                & i_aggregate, i_clw, i_ciw, i_totalice, get_d_dia, &
!                & psd_modified_gamma, psd_marshall_palmer, psd_field_2005, psd_field_2007 
use mod_mie,  only: mue, n0_fac, onethird, gm32gcm3, pi, i_totalice, get_d_dia, &
               & psd_modified_gamma, psd_marshall_palmer, psd_field_2005, psd_field_2007
!INTF_ON

implicit none

!* common variables
real (kind=jprb),    intent ( in)   :: lwc, temp, f_ghz
real (kind=jprb),    intent (inout) :: dia_froz(n_dia)
integer (kind=jpim), intent ( in)   :: psd, dens, particle_habit, i_type, i_scat
logical (kind=jplm), intent ( in)   :: ll_dsd, ll_variable_d
character,           intent ( in)   :: regime
real (kind=jprb),    intent (out)   :: nd(n_dia) 

!INTF_END

!* local variables
real    (kind=jprb):: lwc_int, lwc_fac, n0, d0, lambda, dia_melt, nr, by_dia_density(n_dia)
real    (kind=jprb):: mass(n_dia) 
integer (kind=jpim):: i_dia
real    (kind=jprb):: d_dia_froz
real    (kind=jprb):: a, b, gamma1plusb

#include "predict_psd.interface"
#include "predict_psd_F07.interface"
#include "gamma_dsd.interface"
#include "mod_gamma_dsd.interface"
#include "density_all.interface"

if (psd == psd_field_2005 ) then

  call predict_psd(i_type, lwc, temp, dia_froz, nd, mass, dens, f_ghz, &
    & i_scat, particle_habit, ll_variable_d)

else if (psd == psd_field_2007) then

  call predict_psd_F07(i_type, lwc, temp, f_ghz, dia_froz, nd, mass, dens, regime, &
    & i_scat, particle_habit, ll_variable_d)

else

  call density_all(i_type, i_scat, particle_habit, dens, f_ghz, &
    & dinmeters=dia_froz/100.0_JPRB, density=by_dia_density, a=a, b=b, &
    & gamma1plusb=gamma1plusb)

  if (psd == psd_marshall_palmer) then 

    call gamma_dsd (i_type, temp, lwc, mue (i_type), by_dia_density, n0_fac (i_type), &
      & ll_dsd, n0, d0, lambda, a, b, gamma1plusb)
  
  elseif (psd == psd_modified_gamma) then

    ! Modified gamma (does not yet support variable density)
    if(any(by_dia_density /= by_dia_density(1))) then
      write(*,*) 'Density must be constant with diameter when using modified gamma distributions'
      stop
    endif

  endif

  do i_dia = 1, n_dia

    if (psd == psd_marshall_palmer) then

      ! Marshall Palmer
      nd(i_dia)=n0*(dia_froz(i_dia)**mue(i_type))/exp(lambda*dia_froz(i_dia))

    elseif (psd == psd_modified_gamma) then

      dia_melt = dia_froz (i_dia) * by_dia_density(i_dia) ** onethird
      call mod_gamma_dsd (lwc, dia_melt, nr)
      nd(i_dia) = nr 
 
    endif    

    ! Particle mass [g]
    mass(i_dia) = pi/6.0_jprb*by_dia_density(i_dia)*dia_froz(i_dia)**3.0_jprb
 
  enddo

endif

d_dia_froz = get_d_dia( dia_froz, n_dia )

if ( i_type /= i_totalice ) then

  ! Normalise in case the range of dia_frozen omits part of the size distribution
  lwc_int = sum( nd(:) * mass(:) * d_dia_froz )
  lwc_fac = gm32gcm3 * lwc / lwc_int 
  nd(:)   = lwc_fac * nd(:)

  if ( abs(1.0_JPRB - lwc_fac) > 2.0_JPRB ) then
    write(*,*) 'Renormalisation of greater than 200% suggests a fundamental size distribution problem: ', lwc_fac
    write(*,*) minval(dia_froz), maxval(dia_froz), lwc, i_type
    stop 
  endif 

endif

return
end subroutine set_spectra
