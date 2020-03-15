subroutine mie_sphere (x, m, q_sct, q_ext, g, q_bsct)

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
!           04/09/2008   transferred to fcm (Amy Doherty)
!           09/03/2010   performance enhancements (Alan Geer)
!           28/04/2010   fixed "kind" bug which used to require "-r8" 
!                        compile to prevent convergence failures (Alan Geer)

use parkind1, only: jprb
!INTF_OFF
use parkind1, only: jpim
!INTF_ON
implicit none

!* common variables
real    (kind=jprb), intent ( in) :: x
complex (kind=jprb), intent ( in) :: m
real    (kind=jprb), intent (out) :: q_sct, q_ext, g, q_bsct

!INTF_END      

!* local variables  
real     (kind=jprb):: n1, n2, d0_re, d0_im, f_n, g_n, nn
complex  (kind=jprb):: mc, mx, c_bsct, anbn, bnan
integer  (kind=jpim):: n, n_end

complex (kind=jprb) :: dn, dn_prev, wn, wn_prev1, wn_prev2, an, an_prev, bn, bn_prev

!* how many iterations 
n_end = nint (x + 4.05_jprb * x ** (1.0_jprb / 3.0_jprb) + 2.0_jprb)
if (n_end <    5) n_end =    5
if (n_end > 1000) n_end = 1000

mc = conjg (m)
mx = mc * x
n1 = real (m) * x
n2 = aimag (m) * x
      
d0_re = sin  (n1) * cos  (n1) / (sin (n1) * sin (n1) + sinh (n2) * sinh (n2))
d0_im = sinh (n2) * cosh (n2) / (sin (n1) * sin (n1) + sinh (n2) * sinh (n2))

!* initialize 
dn_prev  = cmplx (d0_re  ,               d0_im, jprb)
wn_prev2 = cmplx (cos (x), -1.0_jprb * sin (x), jprb)
wn_prev1 = cmplx (sin (x),             cos (x), jprb)

q_ext  = 0.0_jprb
q_sct  = 0.0_jprb
q_bsct = 0.0_jprb
g      = 0.0_jprb 
c_bsct = cmplx (0.0_jprb,0.0_jprb, jprb)

do n = 1, n_end + 1

  f_n = 2.0_jprb * n - 1.0_jprb
  g_n = 2.0_jprb * n + 1.0_jprb

  wn = f_n / x * wn_prev1 - wn_prev2
  dn = 1.0_jprb / (n / mx - dn_prev) - n / mx
  
  an = ((dn / mc + n / x) * real (wn) - real (wn_prev1)) / ((dn / mc + n / x) * wn - wn_prev1)
  bn = ((dn * mc + n / x) * real (wn) - real (wn_prev1)) / ((dn * mc + n / x) * wn - wn_prev1)
     
  q_ext  = q_ext  + g_n * (real (an) + real (bn))
  q_sct  = q_sct  + g_n * (abs (an) ** 2.0_jprb + abs (bn) ** 2.0_jprb)
  c_bsct = c_bsct + g_n * (-1.0_jprb) ** n * (an - bn)

  if (n > 1) then
     
    anbn = an_prev * conjg (an) + bn_prev * conjg (bn)
    bnan = an_prev * conjg (bn_prev)
   
    nn = n - 1
    g = g + nn * (nn + 2.0_jprb) / (nn + 1.0_jprb) * real (anbn) + (2.0_jprb * nn + 1.0_jprb) / (nn * (nn + 1.0_jprb)) * real (bnan)

  end if

  dn_prev  = dn
  wn_prev2 = wn_prev1
  wn_prev1 = wn
  an_prev  = an
  bn_prev  = bn

end do

q_ext  = q_ext * 2.0_jprb / (x * x)
q_sct  = q_sct * 2.0_jprb / (x * x)
if (q_sct > q_ext) q_sct = q_ext

g      = g  * 4.0_jprb / (x * x * q_sct)
q_bsct = abs (c_bsct) * abs (c_bsct) / (x * x)

return
end subroutine mie_sphere
