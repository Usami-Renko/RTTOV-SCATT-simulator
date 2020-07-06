subroutine liu_density (liu_habit, a, b, gamma)

! This software was developed within the context of
! the EUMETSAT Satellite Application Facility on
! Numerical Weather Prediction (NWP SAF), under the
! Cooperation Agreement dated 25 November 1998, between
! EUMETSAT and the Met Office, UK, by one or more partners
! within the NWP SAF. The partners in the NWP SAF are
! the Met Office, ECMWF, KNMI and MeteoFrance.
!
! Copyright 2002, EUMETSAT, All Rights Reserved.

! This function implements the ice density relations given in table 1 of 
! Kulie et al. (2010, JAS 67, 3471 - 3487) relating to the Liu (2008, BAMS, 1563-1570)
! ice particle shapes. The habit number is defined in scatdb.c. The return values
! a and b are the coefficients of the mass - particle size relationship:
!
!   m(D) = aD^b
!
! Where D is the maximum dimension of the ice particle and m its mass. Units are SI,
! i.e. m is in kg and D is in m. Kulie et al. assume an ice density of 1000 kg m^-3
! which makes "a" 8% incorrect, but the same numbers are used here for consistency.
!
! For convenience in computing a Gamma size distribution based on these shapes,
! this subtroutine also returns the value of the gamma function for 1+b,
! precomputed offline using IDL's implementation of the gamma function.

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           22/02/2011  initial version (Alan Geer)

use parkind1, only: jprb, jpim

implicit none

! Interface
integer (kind=jpim), intent ( in) :: liu_habit
real    (kind=jprb), intent (out) :: a, b, gamma
!INTF_END

! Data
real (kind=jprb) table_a(0:10), table_b(0:10), gamma_1plusb(0:10)

data table_a(0:10)      / 37.09_jprb, 116.12_jprb, 229.66_jprb, 122.66_jprb, 32.36_jprb, 0.32_jprb,   &
                           0.06_jprb,   0.07_jprb,   0.09_jprb,  0.002_jprb,  0.01_jprb   /
data table_b(0:10)      /  3.00_jprb,   3.00_jprb,   3.00_jprb,   3.00_jprb,  3.00_jprb, 2.37_jprb,   &
                           2.12_jprb,   2.12_jprb,   2.13_jprb,   1.58_jprb,  1.90_jprb   /
data gamma_1plusb(0:10) /  6.00_jprb,   6.00_jprb,   6.00_jprb,   6.00_jprb,  6.00_jprb, 2.8875_jprb, &
                         2.2405_jprb, 2.2405_jprb, 2.2623_jprb, 1.4084_jprb,  1.8274_jprb /

a     = table_a(liu_habit)
b     = table_b(liu_habit)
gamma = gamma_1plusb(liu_habit)

end subroutine liu_density
