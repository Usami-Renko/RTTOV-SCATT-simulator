subroutine mg_ellips (fv_incl, perm_matrix, perm_incl, perm_mix)

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

! Formula for effective permittivity of 2-component media

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)
!           15/03/1013   Added true Maxwell-Garnett formula for reference (Alan Geer)

use parkind1, only: jprb

implicit none

!* Maxwell-Garnett  

!* common variables
real    (kind=jprb), intent ( in) :: fv_incl
complex (kind=jprb), intent ( in) :: perm_matrix, perm_incl
complex (kind=jprb), intent (out) :: perm_mix

!INTF_END

!* local variables
complex (kind=jprb) :: q, gamma

! This is NOT optional, in order to preserve reproducibility of 
! previously-generated Mie Tables.
if (.true.) then

  ! The standard RTTOV-SCATT formula is Fabry & Szyrmer (contrary to earlier 
  ! documentation and code comments!)
  q        = (perm_incl / (perm_incl - perm_matrix)) * log (perm_incl / perm_matrix) - 1.0_jprb
  gamma    = 2.0_jprb * perm_matrix * q / (perm_incl - perm_matrix)

  perm_mix = ((1.0_jprb - fv_incl) * perm_matrix + fv_incl * gamma * perm_incl) &
            /  (1.0_jprb - fv_incl + fv_incl * gamma) 

else

  ! The true Maxwell-Garnett formula. It really doesn't make much difference (Alan Geer)
  q = fv_incl*(perm_incl - perm_matrix)/(perm_incl + 2.0_JPRB*perm_matrix)
  perm_mix = perm_matrix*(1.0_JPRB + 2.0_JPRB*q)/(1.0_JPRB - q)

endif

return
end subroutine mg_ellips
