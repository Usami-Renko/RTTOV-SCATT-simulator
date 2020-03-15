subroutine mie_one_wc ( wavelength, f_ghz, i_type, permwat, ll_melt, dia_froz, &
  & nd, q_ext, q_sct, q_asm, q_bsct, ext, sct, asm, bsct)

! Description:
!
!    Computes bulk scattering parameters for one particular water content, temperature, 
!    frequency and type, given a size distribution, and scattering properties as a function
!    of size.
!
!    See Petty, 2006: "A first course in atmospheric radiation", sect. 12.4
!
!    IN: wavelength -             [cm]
!        f_ghz      - frequency   [GHz]
!        i_type     - hydrometeor type (see mod_mie.F90)
!        ll_melt    - do melting layer of Bauer (2001)
!        dia_froz   - particle diameter or maximum dimension [cm]
!        nd         - size distribution [cm^-4]
!        q_ext      - Extinction cross-section [cm^2]
!        q_sct      - Scattering cross-section [cm^2]
!        q_bsct     - Backscattering cross-section [cm^2]
!        q_asm      - Asymmetry paramter       [ ]
!
!    OUT: ext       - Bulk extinction [cm^-1]
!         sct       - Bulk scattering [cm^-1]
!         bsct      - Bulk backscatter [cm^-1]
!         asm       - Bulk asymmetry   [ ] (weighted by scattering cross-section)
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

! History:
! Version   Date        Comment
! -------   ----        -------
!           10/03/2010  extracted from create_tables_spectra 
!                       to aid parallelisation (Alan Geer)
!           03/03/2011  now just does the integration. Documentation (Alan Geer)

use parkind1, only: jprb, jpim
use mod_mie,  only: n_dia
!INTF_OFF
use mod_mie,  only: i_snow, i_graupel, i_aggregate, i_rain, get_d_dia
!INTF_ON
implicit none

! Interface
real (kind=jprb),    intent (in   ) :: wavelength, f_ghz
integer (kind=jpim), intent (in   ) :: i_type, permwat
logical,             intent (in   ) :: ll_melt
real (kind=jprb),    intent (in   ) :: dia_froz(n_dia), nd(n_dia) 
real (kind=jprb),    intent (in   ) :: q_ext(n_dia), q_sct(n_dia), q_asm(n_dia), q_bsct(n_dia)
real (kind=jprb),    intent (  out) :: ext, sct, asm, bsct
!INTF_END

! Local variables
real    (kind=jprb) :: itgr
integer (kind=jpim) :: i_dia
complex (kind=jprb) :: perm_froz, perm_wat
real    (kind=jprb) :: d_dia_froz

!* Interface blocks
#include "permittivity.interface"
#include "melting_layer.interface"

ext  = 0.0_jprb
sct  = 0.0_jprb
asm  = 0.0_jprb
bsct = 0.0_jprb

if (ll_melt .and. (i_type == i_snow .or. i_type == i_graupel .or. i_type == i_aggregate) ) then

  perm_wat  = permittivity(i_rain, f_ghz, 273.0_jprb, ipermwat=permwat)
  perm_froz = permittivity(i_type, f_ghz, 273.0_jprb)

  call melting_layer (i_type, wavelength, nd, dia_froz, perm_wat, perm_froz, ext, sct, asm, bsct)

else 

  d_dia_froz = get_d_dia( dia_froz, n_dia )

  do i_dia = 1, n_dia       

    ! rest of integrand 
    itgr = nd (i_dia) * d_dia_froz

    ext  = ext  + q_ext(i_dia)                * itgr
    sct  = sct  + q_sct(i_dia)                * itgr
    asm  = asm  + q_asm(i_dia) * q_sct(i_dia) * itgr
    bsct = bsct + q_bsct(i_dia)               * itgr

  end do 

end if

if (sct > ext) sct = ext
if (sct > 0.0_jprb) asm = asm / sct
      
return
end subroutine mie_one_wc
