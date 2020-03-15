subroutine convert_mietable(input_file, ext_tab, ssa_tab, asm_tab)

! Description:
!
!    Creates bulk-scattering coefficient files for use by RTTOV_SCATT
!
!    IN: input_file - input configuration file
!        ext_tab    - Extinction               [km^-1]
!        ssa_tab    - Single scattering albedo [ ]
!        asm_tab    - Asymmetry paramter       [ ]
!
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
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------
!           06/07/2011  Flexible length config section. Add comments (Alan Geer)

use parkind1, only: jprb
!INTF_OFF
use parkind1, only: jpim
use mod_mie, only: n_temp, n_lwc, conv_a, conv_b, &
 & wc_slope, wc_offset, temp_offset, ctype, lmax
!INTF_ON
implicit none

character(len=*), intent(in) :: input_file
real (kind=jprb), dimension (:,:,:,:), intent(in) :: ext_tab, ssa_tab, asm_tab
!INTF_END

integer(kind=jpim), parameter :: nmaxconfig = 100 
character (len=8) :: csat_in, csens_in  
character (len=lmax) :: cconfig(nmaxconfig), cin, cfile  

integer, allocatable :: id_freq(:)
real   , allocatable :: freq_c(:), freq_o(:)
integer (kind=jpim) :: n_freq, n_sensor
integer (kind=jpim) :: m_type, m_freq, i, nconfig
integer (kind=jpim), allocatable, dimension (:) :: id_type 
integer, dimension (8) :: d

character (len=lmax) :: cdir_table

integer :: ifreq, jfreq, itype, itemp, isensor, nside, ii, jj, i_plat, i_inst, iostat

!-------------------------------------------------------------------------------

open(11, file = input_file, status = 'old')

do ! Find header
  read (11,*) cin
  if (trim(cin) == 'HEADER') exit
enddo
read (11,*) cin
read (11,'(a500)') cdir_table
read (11,*) n_freq
read (11,*) n_sensor
read (11,*) m_type 
allocate(id_type(m_type), id_freq(10*n_freq), freq_c(n_freq), freq_o(n_freq))
read (11,*) id_type 
read (11,*) cin

nconfig = 0
do ! Find frequencies and read config lines
  read (11,'(a500)') cin
  if (trim(cin) == 'FREQUENCIES') then
    exit
  else
    nconfig = nconfig+1
    if (nconfig <= nmaxconfig) then
      cconfig(nconfig) = cin
    else
      write(*,*) 'Too many config lines'
      stop
    endif
  endif
enddo
read (11,*) cin

do ifreq = 1, n_freq
   read (11,*) jfreq, freq_c (ifreq), nside, freq_o (ifreq)
end do
do i = 1, 5
   read (11,*) cin
end do

do isensor = 1, n_sensor

  read (11,1110) csat_in, csens_in, i_plat, i_inst
  read (11,*)    m_freq, (id_freq (ii), ii = 1, m_freq)   

  cfile = trim(cdir_table) // '/mietable_' // trim(csat_in) //  &
                      & '_' // trim(csens_in) // '.dat'
  open (20, file = cfile, iostat=iostat)
  if(iostat /= 0) then
    write(*,*)' Problem opening '//trim(cfile)
    exit
  endif

  !* Header
  write (20,'(a16,a8)') ' ! Mie-tables   ', csens_in
  call write_divider(20)

  !* Identification section
  write (20,'(a)') 'IDENTIFICATION'
  write (20,'(a)') ' !'
  write (20,'(2i3,a)') i_plat, i_inst,'             ! platform instrument'
  write (20,'(2a8)') csat_in, csens_in
  write (20,'(a)') ' mw                ! sensor type [ir,mw,hi]' 
  write (20,'(a)') ' 10                ! RTTOV compatibility version'      
  write (20,'(a)') '  2                ! version'
  call date_and_time(VALUES=d)
  write(20, '(i5,i3,i3,a)') d(1), d(2), d(3), '        ! creation date'
  call write_divider(20)

  !* Dimensions section
  write (20,'(a)') 'DIMENSIONS'
  write (20,'(a)') ' !'               
  write (20,'(i3,i2,i3,i4,a)') m_freq, m_type, n_temp, n_lwc, '        ! frequencies hydrometeor-types temperatures water-contents'
  call write_divider(20)

  !* Frequencies
  write (20,'(a)') 'FREQUENCIES'
  write (20,'(30f9.4)') (freq_c (id_freq (jj)), jj = 1, m_freq)
  call write_divider(20)
  
  !* Hydrometeor types  
  write (20,'(a)') 'HYDROMETEOR'    
  write (20,'(100a)') (trim(ctype(id_type(ii)))//' ', ii = 1, m_type)
  call write_divider(20)  
  
  !* Conversions 
  write (20,'(a)') 'CONVERSIONS'  
  write (20,'(a)') ' !  a    b         ! RR = a * LWC^b, [RR]=mm/h, [LWC]=g/m^3'
  do itype = 1, m_type         
    write (20,'(2F5.2,a)') conv_a(id_type(itype)), conv_b(id_type(itype)), '         ! '//ctype(id_type(itype))
  enddo  
  write (20,'(a)') ' ! offset          ! T [K] = index [1..70] + offset  '
  do itype = 1, m_type         
    write (20,'(I3,a)') int(temp_offset(id_type(itype))), '                ! '//ctype(id_type(itype))
  enddo 
  write (20,'(a)') ' ! slope offset    ! LWC [g/m^3] = 10.0 ** (slope * (index [1..nwc] - offset))'
  write (20,'(F4.2,I5)' ) wc_slope, wc_offset             
  call write_divider(20)  

  !* Note the configuration options used in the Mie table calculation (text copied straight from channels.dat) 
  do i=1,nconfig
    write (20,'(a)') trim(cconfig(i))  
  enddo

  !* Write selected output    
  write (20,'(a)') 'EXTINCTION'
  do ifreq = 1, m_freq 
    do itype = 1, m_type
      do itemp = 1, n_temp
        write (20,1111) ext_tab (id_freq (ifreq),itype,itemp,:)
      end do
    end do
  end do

  write (20,'(a)') 'ALBEDO'
  do ifreq = 1, m_freq 
    do itype = 1, m_type
      do itemp = 1, n_temp
        write (20,1111) ssa_tab (id_freq (ifreq),itype,itemp,:)
      end do
    end do
  end do

  write (20,'(a)') 'ASYMMETRY'
  do ifreq = 1, m_freq 
    do itype = 1, m_type
      do itemp = 1, n_temp
        write (20,1111) asm_tab (id_freq (ifreq),itype,itemp,:)
      end do
    end do
  end do

  write (*,*) 'Written: ', trim(cfile)

  close (20)
end do

1110    format (2a8,50i4)
1111    format (5(1x,e23.16))

return

contains

subroutine write_divider(unit)
integer, intent(in) :: unit
write(unit,'(a)')' ! ------------------------------------------------------'
end subroutine write_divider

end subroutine convert_mietable


