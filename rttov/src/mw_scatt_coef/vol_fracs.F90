subroutine vol_fracs (dense_ice, dense_air, dense_snow, fv_wat, fv_ice, fv_air)

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

!* ice/air volume fractions as a function of melt-water volume fraction	

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008   transferred to fcm (Amy Doherty)

use parkind1, only: jprb
     
implicit none
     
!* common variables
real (kind=jprb), intent ( in) :: dense_ice, dense_air, dense_snow, fv_wat
real (kind=jprb), intent (out) :: fv_ice, fv_air

!INTF_END
      
!* local variables
real (kind=jprb) :: fv_ii, fv_aa
      
!* dry
fv_ii = (dense_snow - dense_air) / (dense_ice - dense_air)
fv_aa = 1.0_jprb - fv_ii
            
!* melting
fv_ice = fv_ii * (1.0_jprb - fv_wat)
fv_air = fv_aa * (1.0_jprb - fv_wat)      
      
return
end subroutine vol_fracs
