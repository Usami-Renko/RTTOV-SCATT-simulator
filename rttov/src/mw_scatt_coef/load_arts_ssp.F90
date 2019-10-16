subroutine load_arts_ssp (m_type, i_scat, particle_shape, &
                          ssp_arts, f_arts, T_arts, D_arts)

! Copyright:
!
!    Copyright 2017, Chalmers University of Technology/EUMETSAT, All Rights Reserved.

! History:
! Version   Date        Comment
! -------   ----        -------
!           31/10/2017  initial version (Jana Mendrok)

use parkind1, only: jprb, jpim
use mod_arts, only: nf_max_arts, nT_max_arts, nD_max_arts
!INTF_OFF
use mod_mie, only:  i_arts, lmax
use mod_arts, only: arts_folder, arts_files
!INTF_ON

implicit none

integer (kind=jpim), intent (in) :: m_type
integer (kind=jpim), intent (in) :: i_scat(m_type), particle_shape(m_type)
real (kind=jprb), intent (inout) :: ssp_arts(nf_max_arts, nT_max_arts, nD_max_arts,4,m_type)
real (kind=jprb), intent (inout) :: f_arts(0:nf_max_arts), T_arts(0:nT_max_arts), D_arts(0:nD_max_arts)

!INTF_END

! Local variables
integer (kind=jpim) :: nf, nT, nD, fi, Ti, Di, i_type, ios, fid=16
character(len=lmax) :: line
character(len=*), parameter :: err_base = 'Aborting due to problem in ARTS SSP database: '

!FIXME:
! - check that grids are strictly monotone.


do i_type=1,m_type
  ! only in case ARTS-SSP data is used for this hydormeteor, we read the ARTS-SSP data from an rssp habit file
  if (i_scat(i_type) == i_arts) then
    open(unit=fid, status = 'old', iostat=ios, &
         file = trim(arts_folder)//trim(arts_files(particle_shape(i_type))) )
    if (ios /= 0) then
      write(0,*) err_base
      write(0,*) 'Scattering file not found: '//trim(arts_folder)//trim(arts_files(particle_shape(i_type)))
      stop
    endif
    
    ! parse file
    read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) nf, nT, nD
    if (nf>nf_max_arts) then
      write(0,*) err_base
      write(0,*) 'Number of frequencies of m_type(', i_type, ') (nf=', nf, &
                 ') exceeds allowed maximum number of frequencies', &
                 ' in ARTS-SSP data (nf_max=', nf_max_arts, ').'
      stop
    endif
    f_arts(0) = REAL(nf)
    if (nT>nT_max_arts) then
      write(0,*) err_base
      write(0,*) 'Number of temperature of m_type(', i_type, ') (nT=', nf, &
                 ') exceeds allowed maximum number of temperatures', &
                 ' in ARTS-SSP data (nT_max=', nT_max_arts, ').'
      stop
    endif
    T_arts(0) = REAL(nT)
    if (nD>nD_max_arts) then
      write(0,*) err_base
      write(0,*) 'Number of particle sizes of m_type(', i_type, ') (nD=', nD, &
                 ') exceeds allowed maximum number of particle sizes', &
                 ' in ARTS-SSP data (nD_max=', nD_max_arts, ').'
      stop
    endif
    D_arts(0) = REAL(nD)
    read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) f_arts(1:nf)
    read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) T_arts(1:nT)
    read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) D_arts(1:nD)
    ! for some reason, the read-in of the D causes a linefeed.
    ! so, don't linefeed once more to skip that comment, as we have that skipped implicitly already.
    !read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) line ! Dveq grid. Replace by reading into a REAL array when needed.
    read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) line ! mass grid. Replace by reading into a REAL array when needed.
    read(fid, *, iostat=ios) line ! comment line
    read(fid, *, iostat=ios) line ! m-D's a and b. Replace by reading into REAL scalars when needed.
    read(fid, *, iostat=ios) line ! comment line
    do fi=1,nf
      do Ti=1,nT
        do Di=1,nD
          read(fid, *, iostat=ios) ssp_arts(fi,Ti,Di,:,i_type)
        enddo
      enddo
    enddo
    !fi=1
    !Ti=1
    !do Di=1,nD-1
    !  write(*,*) 'loaded ext(fi=',fi,',Ti=',Ti,'Di=',Di+1,')=',ssp_arts(fi,Ti,Di+1,1,1)
    !enddo
    close(fid)
  end if
enddo


end subroutine load_arts_ssp
