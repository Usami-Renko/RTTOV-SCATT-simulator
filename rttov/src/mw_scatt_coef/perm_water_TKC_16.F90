subroutine perm_water_TKC_16 (f_ghz, t_k, perm_re, perm_im)

! Calculates the permittivity of water following Turner et al. 2016 
! "An improved liquid water absorption model at microwave frequencies for supercooled liquid water clouds"
! Journal of atmospheric and oceanic technology, Volume 33, pages 33-44

! History:
! Version   Date        Comment
! -------   ----        -------
!           13/06/2016  Coded for fortran and rttov (Katrin Lonitz)


use parkind1, only: jprb
!INTF_OFF
use mod_mie,  only: pi
!INTF_ON

implicit none

!* common variables
real (kind=jprb), intent ( in) :: f_ghz, t_k
real (kind=jprb), intent (out) :: perm_re, perm_im

!INTF_END

!* local variables
real (kind=jprb)  :: freq, temp
real (kind=jprb)  :: t_c
real (kind=jprb)  :: a_1, a_2
real (kind=jprb)  :: b_1, b_2
real (kind=jprb)  :: c_1, c_2
real (kind=jprb)  :: d_1, d_2
real (kind=jprb)  :: s_0, s_1, s_2, s_3
real (kind=jprb)  :: eps_s
real (kind=jprb)  :: delta_1, delta_2
real (kind=jprb)  :: tau_1, tau_2
real (kind=jprb)  :: term1_1, term1_2
real (kind=jprb)  :: term2_1, term2_2


!*** The fitting formula works with
!    degr. Celsius, Hz, etc.) so check
!    inputs and convert if neccessary:

!*** Convert frequency to SI
freq = f_ghz * 1e9_jprb 

if (freq > 500e+9_jprb) write(0,*) '!!!Frequency range: 0.5-500 GHz!!!, EXTRAPOLATION'

!*** Convert Temperature from Kelvin to degree Celsius
temp = t_k - 273.15_jprb

!--------------------------------------------------------------------------------------------------------
!COEFFICCIENTS
!--------------------------------------------------------------------------------------------------------
  
  a_1 = 8.111E+1_jprb
  b_1 = 4.434E-3_jprb
  c_1 = 1.302E-13_jprb
  d_1 = 6.627E+2_jprb
  a_2 = 2.025E+0_jprb
  b_2 = 1.073E-2_jprb
  c_2 = 1.012E-14_jprb
  d_2 = 6.089E+2_jprb
  t_c = 1.342E+2_jprb


  s_0 = 8.7914E+1_jprb
  s_1 = -4.0440E-1_jprb
  s_2 = 9.5873E-4_jprb
  s_3 = -1.3280E-6_jprb

!--------------------------------------------------------------------------------------------------------
!CALCULATION
!--------------------------------------------------------------------------------------------------------
 
!Static dielectric permettivity (for all the same)
eps_s = s_0 + s_1 * temp + s_2 * temp**2.0_jprb + s_3 * temp**3.0_jprb

delta_1 = a_1 * EXP(-b_1 * temp)
delta_2 = a_2 * EXP(-b_2 * temp)

tau_1   = c_1 * EXP(d_1 / (temp + t_c))
tau_2   = c_2 * EXP(d_2 / (temp + t_c))

!Relaxation terms
term1_1 = (tau_1**2.0_jprb*delta_1) / (1.0_jprb + (2.0_jprb*pi*freq*tau_1)**2.0_jprb)
term1_2 = (tau_2**2.0_jprb*delta_2) / (1.0_jprb + (2.0_jprb*pi*freq*tau_2)**2.0_jprb) 

term2_1 = (tau_1 * delta_1) / (1.0_jprb + (2.0_jprb*pi*freq*tau_1)**2.0_jprb)
term2_2 = (tau_2 * delta_2) / (1.0_jprb + (2.0_jprb*pi*freq*tau_2)**2.0_jprb)

!permittivity coefficients
perm_re = eps_s - ((2.0_jprb*pi*freq)**2.0_jprb)*(term1_1 + term1_2)
perm_im = (2.0_jprb*pi*freq) * (term2_1 + term2_2)

return
end subroutine perm_water_TKC_16 
