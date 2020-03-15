! Description:
!> @file
!!   Executable for converting an ASCII RTTOV-SCATT Mietable file
!!   to binary format.
!
!> @brief
!!   Executable for converting an ASCII RTTOV-SCATT Mietable file
!!   to binary format.
!!
!! @details
!!   Usage:
!!   $ rttov_ascii2bin_scattcoef.exe \-\-coef-in  ... \-\-coef-out  ...
!!
!!   The output filename should end in '.bin'
!!
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
PROGRAM rttov_ascii2bin_scattcoef

#include "throw.h"
  Use parkind1, Only : jpim, jprb
  Use rttov_types, Only : rttov_coefs, rttov_scatt_coef, rttov_options_scatt
  Use rttov_const, Only : rttov_magic_string, rttov_magic_number
  Use rttov_getoptions
  Use rttov_unix_env, Only : rttov_iargc, rttov_exit

  Implicit None

#include "rttov_opencoeff.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_errorreport.interface"

  ! Local variables
  !-------------------
  ! coefficient structure
  
  Type(rttov_coefs)      :: coef_rttov        ! clear-sky coefficients
  Type(rttov_scatt_coef) :: coef_scatt        ! coefficients

  ! Logical units for input/output
  Integer(Kind=jpim) :: file_id
  ! error return code for subroutines
  Integer(Kind=jpim) :: err,io_status

  Type(rttov_options_scatt) :: opts_scatt
  Character (len=256)   :: f_coef_in = '', f_coef_out = ''

  !- End of header --------------------------------------------------------

  TRY

  If( rttov_iargc() .eq. 0 ) Then
    Print *, "Usage: --coef-in  ... --coef-out  ... "
    Stop
  EndIf

  Call getoption( "--coef-in",  f_coef_in  )
  Call getoption( "--coef-out", f_coef_out )

  file_id = 10
  Call rttov_opencoeff (err, f_coef_in, file_id)
  THROWM(err /= 0, 'RTTOV_ASCII2BIN_SCATTCOEF: Error opening input file '//Trim(f_coef_in))
  INFO('Reading from ASCII coefficient file '//f_coef_in)
  Call rttov_read_scattcoeffs(err, opts_scatt, coef_rttov, coef_scatt, file_id)
  THROWM(err /= 0, 'RTTOV_ASCII2BIN_SCATTCOEF: Error reading file '//Trim(f_coef_in))
  Close ( unit = file_id )

  coef_scatt%conv_rain(:) = 1._JPRB/coef_scatt%conv_rain(:)
  coef_scatt%conv_sp  (:) = 1._JPRB/coef_scatt%conv_sp  (:)
  coef_scatt%scale_water = 1._JPRB/coef_scatt%scale_water
  coef_scatt%offset_water = - coef_scatt%offset_water

  INFO('Writing to binary coefficient file '//f_coef_out)
  Call rttov_opencoeff (err, f_coef_out, file_id, for_output=.True., lbinary=.True.)
  THROWM(err /= 0, 'RTTOV_ASCII2BIN_SCATTCOEF: Error opening output file '//Trim(f_coef_out))

  ! Write a string that could be displayed
  ! Write a real number to be able to check single/double precision
  Write(file_id, iostat=io_status) rttov_magic_string, rttov_magic_number

  Write(file_id, iostat=io_status) &
         & coef_scatt%mfreqm,      &
         & coef_scatt%mtype,       &
         & coef_scatt%mtemp,       &
         & coef_scatt%mwc        

  Write(file_id, iostat=io_status) coef_scatt%mie_freq (:)

  Write(file_id, iostat=io_status)          &
         & coef_scatt%conv_rain(:),         &
         & coef_scatt%conv_sp(:),           &
         & coef_scatt%conv_liq(:),          &
         & coef_scatt%conv_ice(:),          &
         & coef_scatt%conv_totalice(:),     &
         & coef_scatt%offset_temp_rain,     &
         & coef_scatt%offset_temp_sp,       &
         & coef_scatt%offset_temp_liq,      &
         & coef_scatt%offset_temp_ice,      &
         & coef_scatt%offset_temp_totalice, &
         & coef_scatt%scale_water,          &
         & coef_scatt%offset_water

  Write(file_id, iostat=io_status) coef_scatt % ext(:,:,:,:)
  
  Write(file_id, iostat=io_status) coef_scatt % ssa(:,:,:,:)

  Write(file_id, iostat=io_status) coef_scatt % asp(:,:,:,:)

  Close ( unit = file_id )

  THROWM(io_status /= 0, 'RTTOV_ASCII2BIN_SCATTCOEF: Write IO error')

  INFO('RTTOV_ASCII2BIN_SCATTCOEF: Done')

  PCATCH

End Program rttov_ascii2bin_scattcoef
