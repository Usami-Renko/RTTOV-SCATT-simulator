subroutine scattering ( dia_froz, temp, wavelength, f_ghz, i_type, dens, permwat,&
  & i_scat, particle_habit, is_loaded, &
  & ssp_arts, f_arts, T_arts, D_arts, &
  & q_ext, q_sct, q_asm, q_bsct)

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
!    Calculates scattering parameters as a function of particle size, given temperature, 
!    frequency and hydrometeor type. Scattering paramters may either come from the Mie 
!    theory or in the case of frozen precipitation, from the Liu (2008) DDA tables. 
!
!    IN: dia_froz   - particle diameter or maximum dimension [cm]
!        temp       - temperature [K]
!        wavelength -             [cm]
!        f_ghz      - frequency   [GHz]
!        i_type     - hydrometeor type (see mod_mie.F90)   
!        dens       - density parametrization (totalice only)
!        i_scat     - number indicates use of Mie, Liu (2008) DDA or ARTSDB scattering params 
!        permwat    - liquid water permittivity scheme (Liquid Mie spheres only)
!        particle_habit  - habit for above, as defined in Liu (2008) or Baran ens shape 1-6
!        is_loaded  - used by liu stuff
!       
!    OUT: q_ext  - Extinction cross-section [cm^2]
!         q_sct  - Scattering cross-section [cm^2]
!         q_bsct - Backscattering cross-section [cm^2]
!         q_asm  - Asymmetry parameter [ ] 
!

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           02/03/2010  New function (Alan Geer)
!           15/03/2013  Fully-flexible PSD, density & shape options (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)
!           10/01/2018  add option to select liquid water permittivity model (Katrin Lonitz)

use parkind1, only: jprb, jpim
use mod_mie,  only: n_dia
!INTF_OFF
use mod_mie,  only: pi, i_mie, i_liu, i_arts
!INTF_ON
use mod_arts, only: nf_max_arts, nT_max_arts, nD_max_arts

implicit none

! Interface
real (kind=jprb),    intent (in   ) :: temp, wavelength, f_ghz, dia_froz(n_dia)
integer (kind=jpim), intent (in   ) :: i_type, dens, particle_habit
integer (kind=jpim), intent (in   ) :: permwat
integer (kind=jpim), intent (inout) :: is_loaded
integer (kind=jpim), intent (in   ) :: i_scat
real (kind=jprb),    intent (  out) :: q_ext(n_dia), q_sct(n_dia), q_asm(n_dia), q_bsct(n_dia)
real (kind=jprb),    intent (in   ) :: ssp_arts(nf_max_arts, nT_max_arts, nD_max_arts,4)
real (kind=jprb),    intent (in   ) :: f_arts(0:nf_max_arts), T_arts(0:nT_max_arts), D_arts(0:nD_max_arts)
!INTF_END

! Local variables
real    (kind=jprb) :: itgr, x
integer (kind=jpim) :: i_dia
complex (kind=jprb) :: m, perm
real (kind=jprb)    :: density(n_dia)

!* Interface blocks
#include "mie_sphere.interface"
#include "permittivity.interface"
#include "liu_dda.interface"
#include "density_all.interface"
#include "arts_scat.interface"

real (kind=jprb) :: Dinmeters

call density_all(i_type, i_mie, -1_jpim, dens, f_ghz, &
  & dinmeters=dia_froz/100.0_jprb, density=density)

do i_dia = 1, n_dia       

  if (i_scat == i_mie ) then

    ! Normal Mie sphere approach
    perm = permittivity(i_type, f_ghz, temp, ipermwat=permwat, density=density(i_dia)) 
    m = sqrt (perm)
    x = pi * dia_froz(i_dia) / wavelength

    call mie_sphere (x, m, q_sct(i_dia), q_ext(i_dia), q_asm(i_dia), q_bsct(i_dia)) 

    ! Convert efficiencies (dimensionless) to cross sections (cm^2) 
    itgr = pi / 4.0_jprb * dia_froz (i_dia) ** 2.0_jprb 
    q_sct(i_dia)  = q_sct(i_dia) * itgr
    q_ext(i_dia)  = q_ext(i_dia) * itgr
    q_bsct(i_dia) = q_bsct(i_dia) * itgr

  elseif (i_scat == i_liu) then

    ! Use Liu (2008) DDA approximations for snow scattering parameters
    Dinmeters = dia_froz(i_dia)/100.0_jprb
    call liu_dda(f_ghz, temp, particle_habit, Dinmeters, q_ext(i_dia), q_sct(i_dia), q_asm(i_dia), q_bsct(i_dia), is_loaded)

  elseif (i_scat == i_arts) then

    ! Use ARTS SSP database
    Dinmeters = dia_froz(i_dia)/100.0_jprb
    call arts_scat(f_ghz, temp, Dinmeters, &
                   ssp_arts, f_arts, T_arts, D_arts, &
                   q_ext(i_dia), q_sct(i_dia), q_asm(i_dia), q_bsct(i_dia))

  endif

end do

return
end subroutine scattering
