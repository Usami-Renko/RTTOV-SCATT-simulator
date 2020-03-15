subroutine melting_layer (i_type, wavelength, nd, dia_froz, perm_wat, perm_ice, ext, sct, asm, bsct)

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
use mod_mie, only: n_dia
!INTF_OFF
use mod_mie, only: density_wat, density_air, density_ice, density_lmt, density_liq, cd_snow, cd_grpl, &
                   onethird, pi, g, t_ml, t_fl, rh, ka, dv, mw, le, lm, r, va, &
                   i_totalice, i_snow
!INTF_ON
implicit none

!* common variables
real    (kind=jprb), intent ( in), dimension (n_dia) :: dia_froz
real    (kind=jprb), intent ( in), dimension (n_dia) :: nd

integer (kind=jpim), intent ( in) :: i_type
real    (kind=jprb), intent ( in) :: wavelength  
complex (kind=jprb), intent ( in) :: perm_wat, perm_ice
real    (kind=jprb), intent (out) :: ext, sct, asm, bsct

!INTF_END

!* local variables
real (kind=jprb), dimension (n_dia) :: fv_melt, fm_melt, mass_melt

real  (kind=jprb):: es, t_c, t_k, ff, ci, lstar, kappa, nscv, nrel, f_evap, t_fl_c, es_fl, z_km, dz_km, z_m, dz_m, dt_c, x, y
real  (kind=jprb):: density_froz, cd_froz, tvel_froz, mass_froz, nd_froz
real  (kind=jprb):: dia_wat, tvel_wat, cd_wat
real  (kind=jprb):: dia_melt, dia_melt0, d_dia_melt, density_melt, cd_melt, tvel_melt, nd_melt, alpha_melt
real  (kind=jprb):: q_ext, q_sct, q_bsct, q_asm       
real  (kind=jprb):: mult

complex  (kind=jprb):: perm_core, perm_coat, m_core, m_coat

integer (kind=jpim) :: m_lvl, i_lvl, i_dia

#include "perm_melt.interface"
#include "mie_coated_sphere.interface"

!* initialize
ext  = 0.0
sct  = 0.0
asm  = 0.0
bsct = 0.0

!* inital values at FL
fv_melt   (:) = 0.0                                             
fm_melt   (:) = 0.0                                             
mass_melt (:) = 0.0                                             

!* T, es at FL  
t_fl_c = t_fl - 273.15                                     
es_fl  = 6.1078 * exp (17.08085 * t_fl_c / (234.175 + t_fl_c)) 

!* number of levels, z-increment
z_m   = 1000.0
dz_m  = 10.0
z_km  = z_m  * 0.001
dz_km = dz_m * 0.001

m_lvl = nint (z_m / dz_m)
if (m_lvl == 0) m_lvl = 1
dz_m = z_m / m_lvl
dt_c = (t_ml - t_fl) / m_lvl                                     

do i_lvl = 0, m_lvl                                             
  t_k  = t_fl + i_lvl * dt_c                                 
  t_c  = t_k - 273.15                                         
  es   = 6.1078 * exp (17.08085 * t_c / (234.175 + t_c))  

  !* loop over particle size
  dia_melt0  = 0.009
  i_dia      = 0
   
  do i_dia = 1, n_dia
     
    !* frozen precip D, rho, m
    if (i_type == i_totalice .or. i_type == i_snow) then
      density_froz = 0.012 / dia_froz (i_dia)
      cd_froz      = cd_snow
    else        
      density_froz = density_liq (i_type)
      cd_froz      = cd_grpl
    end if
    if (density_froz < density_lmt) density_froz = density_lmt
    if (density_froz > density_ice) density_froz = density_ice
      
    tvel_froz = sqrt (4.0 * g * dia_froz (i_dia) * density_froz / (3.0 * density_air * cd_froz))
    mass_froz = pi / 6.0 * density_froz * dia_froz (i_dia) ** 3.0
   
    !* water D, V, C_dr
    dia_wat  = (dia_froz (i_dia) ** 3.0 * density_froz / density_wat) ** onethird
    tvel_wat = 100.0 * (9.65 - 10.3 / exp (6.0 * dia_wat))
    if (tvel_wat < 100.0) tvel_wat = 100.0
    cd_wat   = 4.0 * g * dia_wat * density_wat / (3.0 * density_air * tvel_wat * tvel_wat)
     
    !* Gamma-DSD
    nd_froz = nd (i_dia)
   
    !* D, rho
    density_melt = density_froz * density_wat / (fm_melt (i_dia) * density_froz + (1.0 - fm_melt (i_dia)) * density_wat)
    dia_melt     = (dia_froz (i_dia) ** 3.0 * density_froz / density_melt) ** onethird
    d_dia_melt   = dia_melt - dia_melt0
      
    !* terminal fall velocity
    cd_melt   = (cd_froz - cd_wat) * (dia_melt - dia_froz (i_dia)) / (dia_froz (i_dia) - dia_wat) + cd_froz
    tvel_melt = sqrt (4.0 * g * dia_melt * density_melt / (3.0 * density_air * cd_melt))

    !* core/coat permittivities   

    call perm_melt (fv_melt (i_dia), fm_melt (i_dia), density_froz, &
                   perm_wat, perm_ice, alpha_melt, perm_core, perm_coat)
    !* refr. indices   
    m_core = sqrt (perm_core)
    m_coat = sqrt (perm_coat)
    !* core/coat size parameters  
    x = pi * dia_melt * alpha_melt / wavelength
    y = pi * dia_melt              / wavelength
    !* Mie            
    call mie_coated_sphere (x, y, m_core, m_coat, q_sct, q_ext, q_asm, q_bsct)
    !* flux conservation      
    nd_melt = nd_froz * tvel_froz / tvel_melt

    mult = pi / 4.0 * nd_melt * dia_melt ** 2.0 * d_dia_melt * dz_km
    ext  = ext  + q_ext         * mult
    sct  = sct  + q_sct         * mult
    asm  = asm  + q_asm * q_sct * mult
    bsct = bsct + q_bsct        * mult

    !* Next melting stage
    !* eqs. 4, 6
    lstar = dia_melt
    nscv  = 0.64
    nrel  = lstar * tvel_melt / va
    kappa = (nscv ** onethird) * sqrt (nrel)

    if (kappa <= 1.0) then
      ff = 1.00 + 0.14 * kappa * kappa
    else
      ff = 0.86 + 0.28 * kappa
    end if
    !* eq. 7
    ci = 0.5 * dia_melt
    !* eq. 2
    f_evap = rh * es * 1000.0 / t_k - es_fl * 1000.0 / t_fl
    if (f_evap < 0.0) f_evap = 0.0

    mass_melt (i_dia) = mass_melt (i_dia) + dz_m / (0.01 * tvel_melt) * 4.0 * pi * ff * ci / lm &
                        * (ka * (t_k - t_fl) + dv * mw * le / r * f_evap)
    if (mass_melt (i_dia) > mass_froz) mass_melt (i_dia) = mass_froz
   
    !* meltwater fractions
    fm_melt (i_dia) = mass_melt (i_dia) / mass_froz
    fv_melt (i_dia) = fm_melt (i_dia) * density_froz / density_wat &
                      / (1.0 - fm_melt (i_dia) * (1.0 - density_froz / density_wat))
      
    dia_melt0 = dia_melt
  end do
end do

!* average over melting layer
ext  = ext  / z_km
sct  = sct  / z_km
asm  = asm  / z_km
bsct = bsct / z_km

return
end subroutine melting_layer
