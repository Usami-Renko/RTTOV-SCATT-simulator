subroutine liu_dda (f_ghz, temp, liu_habit, Dinmeters, q_ext, q_sct, q_asm, q_bsct, is_loaded)

! This software was developed within the context of
! the EUMETSAT Satellite Application Facility on
! Numerical Weather Prediction (NWP SAF), under the
! Cooperation Agreement dated 25 November 1998, between
! EUMETSAT and the Met Office, UK, by one or more partners
! within the NWP SAF. The partners in the NWP SAF are
! the Met Office, ECMWF, KNMI and MeteoFrance.
!
! Copyright 2002, EUMETSAT, All Rights Reserved.

! This is a wrapper to the Liu scatdb function, necessary to deal with:
!
! 1) conversions from Liu units to cm etc. as used through most of the Mie code
! 2) argument kinds are not used by scatdb and the arguments are JPRM, not JPRB floats.
! 
! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           22/02/2011  initial version (Alan Geer)

use parkind1, only: jprb, jpim
!INTF_OFF
use parkind1, only: jprm
!INTF_ON
implicit none

! Interface
integer (kind=jpim), intent ( in) :: liu_habit
integer (kind=jpim), intent (inout) :: is_loaded
real    (kind=jprb), intent ( in) :: f_ghz, temp, Dinmeters        ! [Ghz, K, m] 
real    (kind=jprb), intent (out) :: q_ext, q_sct, q_asm, q_bsct   ! [cm^2]
!INTF_END

real    (kind=jprm) f,t,dmax,cabs,csca,cbsc,g,ph(37),re
integer (kind=jpim) nshp,iret

f    = f_ghz
t    = max( temp, 234.0_JPRB ) ! Liu tables don't go lower than 233.15K.
nshp = liu_habit
dmax = Dinmeters * 1e6_JPRB ! convert to microns

call scatdb(f,t,nshp,dmax,cabs,csca,cbsc,g,ph,re,iret,is_loaded)

if (iret /= 0) then
  write(*,*) 'Aborting due to problem calling Liu(2008) DDA parameterisation: ', iret
  if (f_ghz < 3.0_JPRB) write(*,*) 'This is not available for frequencies below 3 GHz'
  stop
endif

! Convert output [m^2] cross sections to [cm^2] and return extinction
q_ext = (csca + cabs) * 1e4_JPRB
q_sct = csca * 1e4_JPRB
q_bsct = cbsc * 1e4_JPRB

! Output asymmetry parameter is OK (apart from JPRM kind)
q_asm = g

end subroutine liu_dda
