subroutine mie_coated_sphere (x, y, m1, m2, q_sct, q_ext, g, q_bsct)

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

! Calculates scattering and extinction cross sections using the mie
! solution for cloud and precipitating hydrometeors.

! Current Code Owner: SAF NWP

! History:
! Version   Date        Comment
! -------   ----        -------
!           28/04/2010  fixed kind bug which used to require "-r8" compile (Alan Geer)

use parkind1, only: jprb
!INTF_OFF
use parkind1, only: jpim
use mod_mie, only: eps
!INTF_ON
implicit none

!* common variables
complex (kind=jprb), intent ( in) :: m1, m2
real    (kind=jprb), intent ( in) :: x, y
real    (kind=jprb), intent (out) :: q_sct, q_ext, g, q_bsct
      
!INTF_END

!* local variables
real    (kind=jprb) :: psiy_prev2, psiy_prev1, psiy, chiy_prev2, chiy_prev1, chiy
complex (kind=jprb) :: xiy_prev2, xiy_prev1, xiy
complex (kind=jprb) :: chiy2_prev2, chiy2_prev1, chiy2, chix2_prev2, chix2_prev1, chix2
complex (kind=jprb) :: dx1_prev1, dx1, dx2_prev1, dx2, dy2_prev1, dy2
complex (kind=jprb) :: an_prev1, an, bn_prev1, bn 

complex (kind=jprb) :: xm1, xm2, ym1, ym2, m1m2, psix2, psiy2, chipx2, chipy2, chelp1, chelp2
complex (kind=jprb) :: ancap, bncap, dnbar, gnbar, c_bsct

real     (kind=jprb) :: im_xm1, im_xm2, im_ym2, re_chelp1, re_chelp2, rhelp1, rhelp2, dc_bsct
integer  (kind=jpim) :: n_iter, n

xm1  = m1 * x
xm2  = m2 * x
ym1  = m1 * y
ym2  = m2 * y
m1m2 = m2 / m1

!* check convergence
im_xm1 = 0.5_jprb * (xm1 - conjg (xm1))
im_xm2 = 0.5_jprb * (xm2 - conjg (xm2))
im_ym2 = 0.5_jprb * (ym2 - conjg (ym2))
n_iter = nint (y + 4.0_jprb * y ** 0.333_jprb + 2.0_jprb)
if (im_xm1 > 30 .or. im_xm2 > 30 .or. im_ym2 > 30) stop 'coated sphere Mie-routine'

!* initial values
dx1_prev1 = cos (xm1) / sin (xm1)
dx2_prev1 = cos (xm2) / sin (xm2)
dy2_prev1 = cos (ym2) / sin (ym2)

psiy_prev2 =             cos (y)
psiy_prev1 =             sin (y)
chiy_prev2 = -1.0_jprb * sin (y)
chiy_prev1 =             cos (y)
xiy_prev2  = cmplx (psiy_prev2,-chiy_prev2, jprb)
xiy_prev1  = cmplx (psiy_prev1,-chiy_prev1, jprb)

chiy2_prev2 = -1.0_jprb * sin (ym2)
chiy2_prev1 =             cos (ym2)
chix2_prev2 = -1.0_jprb * sin (xm2)
chix2_prev1 =             cos (xm2)

q_ext   = 0.0_jprb
q_sct   = 0.0_jprb
g       = 0.0_jprb
c_bsct  = cmplx (0.0_jprb,0.0_jprb,jprb)
dc_bsct = 0.0_jprb
n       = 0

!* upward recursion
do while (n < n_iter .and. abs (1.0_jprb - dc_bsct) > eps)             
  n = n + 1

  psiy = (2.0_jprb * n - 1.0_jprb) / y * psiy_prev1 - psiy_prev2
  chiy = (2.0_jprb * n - 1.0_jprb) / y * chiy_prev1 - chiy_prev2
  xiy  = cmplx (psiy,-chiy, jprb)

  dy2  = 1.0_jprb / (n / ym2 - dy2_prev1) - n / ym2
  dx1  = 1.0_jprb / (n / xm1 - dx1_prev1) - n / xm1
  dx2  = 1.0_jprb / (n / xm2 - dx2_prev1) - n / xm2

  chix2 = (2.0_jprb * n - 1.0_jprb) / xm2 * chix2_prev1 - chix2_prev2
  chiy2 = (2.0_jprb * n - 1.0_jprb) / ym2 * chiy2_prev1 - chiy2_prev2
     
  chipx2 = chix2_prev1 - n * chix2 / xm2
  chipy2 = chiy2_prev1 - n * chiy2 / ym2

  psix2 = 1.0_jprb / (chix2 * dx2 - chipx2)  
  psiy2 = 1.0_jprb / (chiy2 * dy2 - chipy2)

  ancap = psix2 * (m1m2 * dx1 - dx2 ) / (m1m2 * dx1 * chix2 - chipx2)        
  bncap = psix2 * (m1m2 * dx2 - dx1 ) / (m1m2 * chipx2 - chix2 * dx1 )

  dnbar = (dy2 - ancap * chipy2 / psiy2) / (1.0_jprb - ancap * chiy2 / psiy2)
  gnbar = (dy2 - bncap * chipy2 / psiy2) / (1.0_jprb - bncap * chiy2 / psiy2)   
   
  an = ((dnbar / m2 + n / y) * psiy - psiy_prev1) &
     / ((dnbar / m2 + n / y) * xiy  - xiy_prev1)
      
  bn = ((gnbar * m2 + n / y) * psiy - psiy_prev1) &
     / ((gnbar * m2 + n / y) * xiy  - xiy_prev1)

  q_sct  = q_sct  + (2.0_jprb * n + 1.0_jprb) * (abs (an) * abs (an) &
                  +  abs (bn) * abs (bn))
  q_ext  = q_ext  + (2.0_jprb * n + 1.0_jprb) * (real (an) + real (bn))
  c_bsct = c_bsct + (2.0_jprb * n + 1.0_jprb) * (-1.0_jprb) ** n * (an - bn)     

  if (n > 1) then   
    chelp1 = an_prev1 * conjg (an) + bn_prev1 * conjg (bn)
    chelp2 = an_prev1 * conjg (bn_prev1)

    rhelp1  = (n * n - 1.0_jprb) / n 
    rhelp2  = (2.0_jprb * n - 1.0_jprb) / (n * n - n)
     
    re_chelp1 = 0.5_jprb * (chelp1 + conjg (chelp1))
    re_chelp2 = 0.5_jprb * (chelp2 + conjg (chelp2))

    g = g + rhelp1 * re_chelp1 + rhelp2 * re_chelp2
  end if

  dc_bsct = abs ((c_bsct - (2.0_jprb * n + 1.0_jprb) * (-1.0_jprb) ** n * (an - bn)) / c_bsct)

  dx1_prev1 = dx1
  dx2_prev1 = dx2
  dy2_prev1 = dy2

  psiy_prev2 = psiy_prev1
  psiy_prev1 = psiy
  chiy_prev2 = chiy_prev1
  chiy_prev1 = chiy
  xiy_prev2  = xiy_prev1
  xiy_prev1  = xiy

  chiy2_prev2 = chiy2_prev1
  chiy2_prev1 = chiy2
  chix2_prev2 = chix2_prev1
  chix2_prev1 = chix2

  an_prev1 = an
  bn_prev1 = bn 

end do

q_ext  = q_ext * 2.0_jprb / (y * y)
q_sct  = q_sct * 2.0_jprb / (y * y)
g      = g * 4.0_jprb / (y * y * q_sct)
q_bsct = abs (c_bsct) * abs (c_bsct) / (y * y)    

return
end subroutine mie_coated_sphere
