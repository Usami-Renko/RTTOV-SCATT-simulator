subroutine arts_scat (f_ghz, temp, Dinmeters, &
                      ssp_arts, f_arts, T_arts, D_arts, &
                      q_ext, q_sct, q_asm, q_bsct)

! Copyright:
!
!    Copyright 2017, Chalmers University of Technology/EUMETSAT, All Rights Reserved.

! All interpolations done here are linear.

! History:
! Version   Date        Comment
! -------   ----        -------
!           31/10/2017  initial version (Jana Mendrok)

use parkind1, only: jprb
use mod_arts, only: nf_max_arts, nT_max_arts, nD_max_arts
!INTF_OFF
use parkind1, only: jpim
use mod_arts, only: f_extpol, T_extpol, D_extpol
!INTF_ON

implicit none

real (kind=jprb), intent (in   ) :: f_ghz, temp, Dinmeters        ! [GHz, K, m] 
real (kind=jprb), intent (in   ) :: ssp_arts(nf_max_arts, nT_max_arts, nD_max_arts,4)
real (kind=jprb), intent (in   ) :: f_arts(0:nf_max_arts), T_arts(0:nT_max_arts), D_arts(0:nD_max_arts) ! [Hz, K, m]
real (kind=jprb), intent (  out) :: q_ext, q_sct, q_asm, q_bsct   ! [cm^2]

!INTF_END

! Local variables
real    (kind=jprb) :: f_hz, fw, Tw, Dw, fw1, Tw1, Dw1, rtmp
integer (kind=jpim) :: nf, nT, nD, fi, Ti, Di
character(len=*), parameter :: err_base = 'Aborting due to problem in ARTS SSP database: '

f_hz = f_ghz*1e9_JPRB ! [Hz]

nf = INT(f_arts(0))
nT = INT(T_arts(0))
nD = INT(D_arts(0))

! JM: 
! It seems more appropriate to do the f&T finding outside the D-loop (maybe even
! up in those loops already. i could re-define the rssp format (and rewrite the
! rssp reading routine) such that passed sub-arrays are contigious...). unless
! a similar loading/caching proedure as for liu/baran is implemented (not by me,
! though.)

! find freq interpolation bounds. allow some extrapolation.
rtmp = (1.+f_extpol)*f_arts(1) - f_extpol*f_arts(2)
if (f_hz<rtmp) then
  write(*,*) err_base
  write(*,'(2(A,F5.1,A))') &
    'Lowest frequency given for current habit is ', &
    f_arts(1)/1e9_JPRB, 'GHz.\n', &
    'Extrapolation set to be allowed at maximum down to ', &
    rtmp/1e9_JPRB, 'GHz.'
  stop
endif

rtmp = (1.+f_extpol)*f_arts(nf) - f_extpol*f_arts(nf-1)
if (f_hz>rtmp) then
  write (*,*) 'hit upper freq bound.'
  write(*,*) err_base
  write(*,'(2(A,F5.1,A))') &
    'Highest frequency given for current habit is ', &
    f_arts(nf)/1e9_JPRB, 'GHz.\n', &
    'Extrapolation set to be allowed at maximum up to ', &
    rtmp/1e9_JPRB, 'GHz.'
  stop
endif

fi = 1
do while ( (fi<nf-1) .and. (f_arts(fi+1)<f_hz) )
  fi = fi+1
enddo


! find temp interpolation bounds. allow some extrapolation.
rtmp = (1.+T_extpol)*T_arts(1) - T_extpol*T_arts(2)
if (temp<rtmp) then
  write(*,*) err_base
  write(*,'(2(A,F5.1,A))') &
    'Lowest temperature given for current habit is ', &
    T_arts(1), 'K.\n', &
    'Extrapolation set to be allowed at maximum down to ', &
    rtmp, 'K.'
  stop
endif

rtmp = (1.+T_extpol)*T_arts(nT) - T_extpol*T_arts(nT-1)
if (temp>rtmp) then
  write(*,*) err_base
  write(*,'(2(A,F5.1,A))') &
    'Highest temperature given for current habit is ', &
    T_arts(nT), 'K.\n', &
    'Extrapolation set to be allowed at maximum up to ', &
    rtmp, 'K.'
  stop
endif

Ti = 1
do while ( (Ti<nT-1) .and. (T_arts(Ti+1)<temp) )
  Ti = Ti+1
enddo


! find size interpolation bounds
rtmp = (1.+D_extpol)*D_arts(1) - D_extpol*D_arts(2)
if (Dinmeters<rtmp) then
  write(*,*) err_base
  write(*,'(2(A,F5.1,A))') &
    'Smallest maximum diameters given for current habit is ', &
    D_arts(1)*1e2_JPRB, 'cm.\n', &
    'Extrapolation set to be allowed at maximum down to ', &
    rtmp*1e2_JPRB, 'cm.'
  stop
endif

rtmp = (1.+D_extpol)*D_arts(nD) - D_extpol*D_arts(nD-1)
if (Dinmeters>rtmp) then
  write(*,*) err_base
  write(*,'(2(A,F5.1,A))') &
    'Largest maximum diameters given for current habit is ', &
    D_arts(nD)*1e2_JPRB, 'cm.\n', &
    'Extrapolation set to be allowed at maximum up to ', &
    rtmp*1e2_JPRB, 'cm.'
  stop
endif

Di = 1
do while ( (Di<nD-1) .and. (D_arts(Di+1)<Dinmeters) )
  Di = Di+1
enddo


! Determine weights
fw = (f_arts(fi+1)-f_hz)      / (f_arts(fi+1)-f_arts(fi))
Tw = (T_arts(Ti+1)-temp)      / (T_arts(Ti+1)-T_arts(Ti))
Dw = (D_arts(Di+1)-Dinmeters) / (D_arts(Di+1)-D_arts(Di))
fw1 = 1.-fw
Tw1 = 1.-Tw
Dw1 = 1.-Dw

! Do the interpolation
q_ext = (ssp_arts(fi,  Ti,  Di,  1)*fw *Tw *Dw  + &
         ssp_arts(fi,  Ti,  Di+1,1)*fw *Tw *Dw1 + &
         ssp_arts(fi,  Ti+1,Di,  1)*fw *Tw1*Dw  + &
         ssp_arts(fi,  Ti+1,Di+1,1)*fw *Tw1*Dw1 + &
         ssp_arts(fi+1,Ti,  Di,  1)*fw1*Tw *Dw  + &
         ssp_arts(fi+1,Ti,  Di+1,1)*fw1*Tw *Dw1 + &
         ssp_arts(fi+1,Ti+1,Di,  1)*fw1*Tw1*Dw  + &
         ssp_arts(fi+1,Ti+1,Di+1,1)*fw1*Tw1*Dw1) *1e4_JPRB ! [m2] -> [cm2] conversion
q_sct = (ssp_arts(fi,  Ti,  Di,  2)*fw *Tw *Dw  + &
         ssp_arts(fi,  Ti,  Di+1,2)*fw *Tw *Dw1 + &
         ssp_arts(fi,  Ti+1,Di,  2)*fw *Tw1*Dw  + &
         ssp_arts(fi,  Ti+1,Di+1,2)*fw *Tw1*Dw1 + &
         ssp_arts(fi+1,Ti,  Di,  2)*fw1*Tw *Dw  + &
         ssp_arts(fi+1,Ti,  Di+1,2)*fw1*Tw *Dw1 + &
         ssp_arts(fi+1,Ti+1,Di,  2)*fw1*Tw1*Dw  + &
         ssp_arts(fi+1,Ti+1,Di+1,2)*fw1*Tw1*Dw1) *1e4_JPRB ! [m2] -> [cm2] conversion
q_asm =  ssp_arts(fi,  Ti,  Di,  3)*fw *Tw *Dw  + &
         ssp_arts(fi,  Ti,  Di+1,3)*fw *Tw *Dw1 + &
         ssp_arts(fi,  Ti+1,Di,  3)*fw *Tw1*Dw  + &
         ssp_arts(fi,  Ti+1,Di+1,3)*fw *Tw1*Dw1 + &
         ssp_arts(fi+1,Ti,  Di,  3)*fw1*Tw *Dw  + &
         ssp_arts(fi+1,Ti,  Di+1,3)*fw1*Tw *Dw1 + &
         ssp_arts(fi+1,Ti+1,Di,  3)*fw1*Tw1*Dw  + &
         ssp_arts(fi+1,Ti+1,Di+1,3)*fw1*Tw1*Dw1
q_bsct= (ssp_arts(fi,  Ti,  Di,  4)*fw *Tw *Dw  + &
         ssp_arts(fi,  Ti,  Di+1,4)*fw *Tw *Dw1 + &
         ssp_arts(fi,  Ti+1,Di,  4)*fw *Tw1*Dw  + &
         ssp_arts(fi,  Ti+1,Di+1,4)*fw *Tw1*Dw1 + &
         ssp_arts(fi+1,Ti,  Di,  4)*fw1*Tw *Dw  + &
         ssp_arts(fi+1,Ti,  Di+1,4)*fw1*Tw *Dw1 + &
         ssp_arts(fi+1,Ti+1,Di,  4)*fw1*Tw1*Dw  + &
         ssp_arts(fi+1,Ti+1,Di+1,4)*fw1*Tw1*Dw1) *1e4_JPRB ! [m2] -> [cm2] conversion

end subroutine arts_scat
