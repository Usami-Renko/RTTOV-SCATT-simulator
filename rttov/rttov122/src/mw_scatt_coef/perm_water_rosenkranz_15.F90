subroutine perm_water_rosenkranz_15 (f_ghz, t_k, perm_re, perm_im)

! Calculates the permittivity of water
! Original source: Rosenkranz "A Model for the Complex Dielectric Constant of
! Supercooled Liquid Water at Microwave Frequencies" IEEE TRANSACTIONS ON
! GEOSCIENCE AND REMOTE SENSING, VOL. 53, NO. 3, MARCH 2015

! History:
! Version   Date        Comment
! -------   ----        -------
!           18/08/2014  Adapted for rttov (Katrin Lonitz)

use parkind1, only: jprb

implicit none

!* common variables
real (kind=jprb), intent ( in) :: f_ghz, t_k
real (kind=jprb), intent (out) :: perm_re, perm_im

!INTF_END

!* local variables
real (kind=jprb)  :: T_C, THETA
real (kind=jprb)  :: EPS_S
real (kind=jprb)  :: DELTA_R, DELTA_B,GAMMA_R, V1
complex (kind=jprb) :: Z2, Z1, Z2_C, Z1_C, F_R, F_B, EPS

!*** The fitting formula works with
!    degr. Celsius and GHz so check
!    inputs and convert if neccessary:

!*** Convert Temperature from Kelvin to degree Celsius
T_C = T_K - 273.15_jprb 

!--------------------------------------------------------------------------------------------------------
!CALCULATION OF eps(f_ghz, t_c) according to eq (2b)
!--------------------------------------------------------------------------------------------------------

!*** Calculate static dielectric constant Eq.(8)
THETA = 300._jprb / T_K
EPS_S = -43.7527_jprb*(THETA**0.05_jprb) + 299.504_jprb*(THETA**1.47_jprb) - &
        399.364_jprb*(THETA**2.11_jprb) + 221.327_jprb*(THETA**2.31_jprb)

!*** Calculate DELTA_R Eq(9)
DELTA_R = 80.69715_jprb*EXP( -T_C / 226.45_jprb)

!*** Calculate DELTA_B Eq(11)
DELTA_B = 4.008724_jprb*EXP( -T_C / 103.05_jprb)

!*** Calculate F_R Eq(6) using Eq(10)
GAMMA_R = 1164.023_jprb*EXP(-651.4728_jprb / (T_C + 133.07_jprb) )
F_R = GAMMA_R / CMPLX(GAMMA_R, F_GHZ, jprb)

!*** Calculate F_B Eq(7) using Eq.(12 -14)
V1 = 10.46012_jprb + 0.1454962_jprb*T_C + 0.063267156_jprb*(T_C**2._jprb) + 0.00093786645_jprb*(T_C**3._jprb)
Z1 = CMPLX(-0.75_jprb*V1,V1, jprb)
Z1_C = CONJG(Z1)

Z2 = CMPLX(-4500._jprb,2000._jprb, jprb)
Z2_C = CONJG(Z2)

! old
!F_B = (LOG( CMPLX(Z2,-F_GHZ, jprb) / CMPLX(Z1,-F_GHZ, jprb)) ) / (2._jprb*LOG(Z2/Z1) ) + &
!    & (LOG( CMPLX(Z2_C,-F_GHZ, jprb) / CMPLX(Z1_C,-F_GHZ, jprb)) ) / (2._jprb*LOG(Z2_C/Z1_C))

!new 
F_B = (LOG( (Z2 - CMPLX(0,F_GHZ, jprb)) / (Z1 - CMPLX(0,F_GHZ, jprb))) ) / (2._jprb*LOG(Z2/Z1) ) + &
    & (LOG( (Z2_C - CMPLX(0,F_GHZ, jprb)) / (Z1_C - CMPLX(0,F_GHZ, jprb))) ) / (2._jprb*LOG(Z2_C/Z1_C))

!*** Calculate permittivity(complex dielectic constant) EQ (2b)
EPS = EPS_S - ( (DELTA_R*(1._jprb - F_R)) + (DELTA_B*(1._jprb - F_B)) )

perm_re = REAL(EPS)
perm_im = -AIMAG(EPS)

return
end subroutine perm_water_rosenkranz_15
