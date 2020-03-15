subroutine predict_psd(i_type,iwc,tk,Dcm,nd,mass,dens,f_ghz,i_scat,particle_habit,ll_variable_d)

! Description:
!   Subroutine copied from Paul fields pvwave code compute_psd.pro
!   This routine uses the second moment and temperature to predict the
!   ice PSD.
!   The moment relations use the 2nd moment as a reference moment. If
!   IWC is computed using a different moment then the second moment 
!   needs to be estimated first.
!   (To avoid this step the moment relations could be found for a different
!   reference moment - requires reanalysis of the original data)
!   See 'Parametrization of ice particle size
!   distributions for mid-latitude stratiform cloud' by Field et al.
!   for more details (submitted to QJRMS August 2004)
!
! Copyright:
!   This software was developed within the context of
!   the EUMETSAT Satellite Application Facility on
!   Numerical Weather Prediction (NWP SAF), under the
!   Cooperation Agreement dated 25 November 1998, between
!   EUMETSAT and the Met Office, UK, by one or more partners
!   within the NWP SAF. The partners in the NWP SAF are
!   the Met Office, ECMWF, KNMI and MeteoFrance.
!
!   Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           04/09/2008  transferred to fcm (Amy Doherty)
!           02/03/2011  normalise (Alan Geer)
!           12/05/2011  Corrected density parameters. Andrew Smith
!           13/03/2013  Fully-flexible PSD, density and shape (Alan Geer)
!           31/10/2017  Adapted for ARTS-SSDB use (Jana Mendrok)
!
! Code Description: 
!   Language:           Fortran 90. 
!   Software Standards: "European Standards for Writing and 
!     Documenting Exchangeable Fortran 90 Code".

use parkind1, only: jprb, jpim, jplm
use mod_mie,  only: n_dia
!INTF_OFF
use mod_mie,  only: not_available
!INTF_ON

implicit none

integer (kind=jpim), intent(in) :: i_type ! Hydrometeor type
real (kind=jprb),    intent(in) :: iwc   ! IWC (g/m^3)
real (kind=jprb),    intent(in) :: tk    ! Temperature (K)
integer (kind=jpim), intent(in) :: dens  ! Density id (1 or 3)

real (kind=jprb), intent(inout) :: Dcm(n_dia)  ! Particle diameter (cm)
real (kind=jprb), intent(out) :: nd(n_dia)     ! Size distribution
real (kind=jprb), intent(out) :: mass(n_dia)

integer  (kind=jpim), intent(in) :: particle_habit
integer  (kind=jpim), intent(in) :: i_scat
real    (kind=jprb), intent ( in) :: f_ghz  ! Frequency, [GHz]
logical (kind=jplm), intent ( in) :: ll_variable_d

!INTF_END

! local variables

! real    (kind=jprb) :: D (n_dia)
real    (kind=jprb) :: dN_dD (n_dia)
real    (kind=jprb) :: xx (n_dia)
real    (kind=jprb) :: phi (n_dia)
real    (kind=jprb) :: tc, x, y, iwcsi, loga_, M3
real    (kind=jprb) :: My, M2, a_, b_
real    (kind=jprb) :: a1, a2, a3, a4, a5, a6, a7, a8, a9, a10
real    (kind=jprb) :: b1, b2, b3, b4, b5, b6, b7, b8, b9, b10
integer (kind=jpim) :: n, i

#include "density_all.interface"

!-------------------------------------------------------------------------------

! Set refactor and exponent x, y in mass-size relation m=xD^y
call density_all(i_type, i_scat, particle_habit, dens, f_ghz, a=x, b=y)

if (x == not_available) then
  write(*,*) 'Invalid density treatment selected for use with Field (2005) PSD'
  stop
endif

! define constants for moment relations based on second moment as reference
!----------------------------------------------------------
a1=5.065339_jprb
a2=-0.062659_jprb
a3=-3.032362_jprb
a4=0.029469_jprb
a5=-0.000285_jprb
a6=0.31255_jprb
a7=0.000204_jprb
a8=0.003199_jprb
a9=0.0_jprb
a10=-0.015952_jprb

b1=0.476221_jprb
b2=-0.015896_jprb
b3=0.165977_jprb
b4=0.007468_jprb
b5=-0.000141_jprb
b6=0.060366_jprb
b7=0.000079_jprb
b8=0.000594_jprb
b9=0.0_jprb
b10=-0.003577_jprb

! check in cloud temperature
!----------------------------------
! first change from Kelvin to degress C

tc = tk - 273.0

if (tc.gt.0) then
  write(*,*) 'in-cloud temperature was too warm (>0C)', tc, tk
  stop
endif
!if (tc.lt.-60) then
!  write(*,*) 'Warning: in-cloud temperature colder than included in original data'
!endif

!------------------------------------

! Calculate second moment
iwcsi = iwc/1E3_jprb !kg/m3
My = iwcsi/x      !yth moment

! Use moment relations to get the other moment required for predicting PSD
M2 = My         ! Default for UM y=2 (do not need to do inversion in 
                !following if-block to avoid slight error introduced by plane fitting)

if (y.ne.2._jprb) then
  n = y
  loga_ = a1+a2*tc+a3*n+a4*tc*n+a5*tc**2+a6*n**2+a7*tc**2*n+a8*tc*n**2+a9*tc**3+a10*n**3
  a_ = 10**loga_
  b_ = b1+b2*tc+b3*n+b4*tc*n+b5*tc**2+b6*n**2+b7*tc**2*n+b8*tc*n**2+b9*tc**3+b10*n**3
  M2 = (My/a_)**(1.0_jprb/b_)
endif
  
! check M2 range (1E-5 to 3E-2 m-1)

if (M2.lt.0.0) then
  write(*,*) 'Negative M2'
  stop
endif

!if (m2.lt.1E-5) then
!  write(*,*)'Warning: M2 outside of range encountered in original data - too small', M2
!endif
!if (m2.gt.3e-2) then
!  write(*,*)'Warning: M2 outside of range encountered in original data - too large', M2
!endif

! Use 2nd and third moment to predict PSD
n=3
loga_=a1+a2*tc+a3*n+a4*tc*n+a5*tc**2+a6*n**2+a7*tc**2*n+a8*tc*n**2+a9*tc**3+a10*n**3
a_=10**loga_
b_=b1+b2*tc+b3*n+b4*tc*n+b5*tc**2+b6*n**2+b7*tc**2*n+b8*tc*n**2+b9*tc**3+b10*n**3
M3=a_*M2**b_

!define dimensionless sizes
do i = 1,n_dia
  if (.not.ll_variable_d) then

    ! Work with existing size ranges
    xx(i) = (m2/m3) * Dcm(i) / 1e2_jprb

  else

    ! Variable size ranges used for totalice
    xx(i) = i/8.0_jprb

  endif
enddo

do i = 1,n_dia

  phi(i)=490.6_jprb*exp(-20.78_jprb*xx(i))+17.46_jprb*xx(i)**0.6357_jprb*exp(-3.29_jprb*xx(i)) !universal distribution

  !;;scale universal psd;;;;;
  dN_dD(i)=phi(i)*m2**4/m3**3 !psd m^-4

  ! Overwrite size ranges 
  Dcm(i)=1e2_jprb*xx(i)*(m3/m2) !particle sizes in cm to match rest of code

enddo

! Mass distribution [g] (x and y are for SI units)
do i = 1, n_dia
  mass(i) = 1e3_jprb * x * (Dcm(i)/100.0_jprb) ** y
enddo

nd = dN_dD*1E-8_jprb ! convert from m^-4 to cm^-4

end subroutine predict_psd
