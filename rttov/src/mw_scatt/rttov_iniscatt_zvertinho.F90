!
Subroutine rttov_iniscatt_zvertinho (&
      & errorstatus,       &! out
      & lradiance,         &! in
      & opts,              &! in
      & nlevels,           &! in
      & nchannels,         &! in
      & nprofiles,         &! in
      & chanprof,          &! in
      & frequencies,       &! in
      & profiles,          &! in  
      & cld_profiles,      &! in 
      & coef_rttov,        &! in
      & nmietables,        &! in
      & coef_scatt_sets,   &! in
      & lshapetable,       &! in
      & transmission,      &! in
      & calcemiss,         &! in
      & usercfrac,         &! in
      & angles,            &! out
      & scatt_aux,         &! inout
      & gh)                 ! in 

  !
  ! Description:
  ! Calculates some variables related to the input precipitation profile
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more painiscattrtners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! 
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002   Initial version     (F. Chevallier)
  !   1.1    05/2003   RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004   Added polarimetry   (R. Saunders)
  !   1.3    08/2004   Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004   Clean-up            (P. Bauer)
  !   1.5    10/2005   Fixes for rttov8 indexing   (U. O'Keeffe)
  !   1.6    11/2005   Add errorstatus to arguments (J. Cameron)   
  !   1.7    11/2007   RTTOV9 / cleanup    (A. Geer)   
  !   1.8    03/2008   Revised cloud partitioning (A. Geer)  
  !   1.9    03/2009   Safety check on cloud fraction (A. Geer) 
  !   1.10   11/2009   User may supply the average cloud fraction (A. Geer)
  !   1.11   11/2009   RTTOV transmittances / optical depths come on levels (A Geer)
  !   1.12   11/2012   Remove "old cloud" option; new subroutine for hydrometeors (A Geer)
  !   1.13   11/2017   R/T now done with radiances, not Tb (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  ! 
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :           &
       & rttov_coef                ,&
       & rttov_scatt_coef          ,&
       & rttov_transmission        ,&
       & rttov_geometry            ,&
       & rttov_profile_scatt_aux   ,&
       & rttov_profile             ,&
       & rttov_profile_cloud       ,&
       & rttov_chanprof            ,&
       & rttov_options

  Use parkind1, Only : jpim, jplm, jprb
!INTF_OFF
  Use rttov_types, Only : rttov_transmission_aux

  Use rttov_const, Only:      &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & gravity,             &
       & pressure_top,        &
       & rgc,                 &
       & mair,                &
       & tcosmic,             &
       & max_scatt_optical_depth

  Use parkind1, Only : jprb, jplm

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Logical (Kind=jplm), Intent (in) :: lradiance ! Computation in radiance, not TB?  
  Type(rttov_options), Intent (in) :: opts      ! RTTOV options
  Integer (Kind=jpim), Intent (in) :: nlevels   ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent(out) :: errorstatus             ! Error return code
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Logical (Kind=jplm), Intent (in) :: calcemiss   (nchannels) ! Emissivity flags
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction
  Type (rttov_profile),             Intent (in)    :: profiles (nprofiles)        ! Atmospheric profiles
  Type (rttov_coef),                Intent (in)    :: coef_rttov                  ! RTTOV Coefficients

  Integer (Kind=jpim),              Intent (in)    :: nmietables                  ! number of mietables
  Type (rttov_scatt_coef),          Intent (in)    :: coef_scatt_sets(nmietables) ! RTTOV_SCATT Coefficients
  Integer (Kind=jpim),              Intent (in)    :: lshapetable(nprofiles, nlevels) ! shape index

  Type (rttov_profile_cloud),       Intent (in)    :: cld_profiles (nprofiles)    ! Cloud profiles
  Type (rttov_transmission),        Intent (in)    :: transmission                ! Transmittances and optical depths
  Type (rttov_geometry),            Intent (out)   :: angles (nprofiles)          ! Zenith angles
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux                   ! Auxiliary profile variables
  Real (Kind=jprb),                 Intent (in)    :: gh(nprofiles, nlevels)      ! Xie: Geopotential Height [km]

!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, ichan, ichanid
  Real    (Kind=jprb) :: dp2dz
  Real    (Kind=jprb) :: z1, z2                                   ! Xie: Geopotential Height [km]

  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf      ! Pressure levels [hPa]
  Real (Kind=jprb), Dimension (nchannels,nlevels)   :: od         ! Single layer optical depths on our levels
  Real (Kind=jprb), Dimension (nlevels+2)           :: od_rttov   ! Single layer optical depths on RTTOV levels
  Real (Kind=jprb), Dimension (nchannels)           :: zod_up_cld ! Optical depth from top of the atmosphere 
  Real (Kind=jprb), Dimension (nprofiles)           :: tbd_prof   ! Boundary temperature, by profile
  Real (Kind=jprb), Dimension (nprofiles,nlevels+1) :: ghfh       ! Xie: Geopotential Height on ph(omit the cosmos level) 
  
  Type (rttov_transmission_aux) :: transmissioncld  ! Clear+cloud transmittances with cloud

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_iniscatt_zvertinho '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
#include "rttov_mieproc_vertinho.interface"
#include "rttov_iniedd.interface"
#include "rttov_hydro.interface"
#include "rttov_calcemis_mw.interface"
#include "rttov_setgeometry.interface"
#include "rttov_errorreport.interface"

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INISCATT_ZVERTINHO',0_jpim,zhook_handle)

  ! write(*,*) "[in]: rttov_iniscatt_vertinho"

  errorstatus = errorstatus_success

  ! NB this only truly needed if using FASTEM 3 
  allocate (transmissioncld % thermal_path1)
  allocate (transmissioncld % thermal_path1 % tau_surf (0:0,nchannels))

  dp2dz = - rgc / gravity / mair

  scatt_aux % ext (:,:) = 0.0_JPRB
  scatt_aux % ssa (:,:) = 0.0_JPRB
  scatt_aux % asm (:,:) = 0.0_JPRB

  ! Get ghfh

  Do iprof = 1, nprofiles

    if (iprof .eq. 1) then
      PRINT *, "GH", gh(1, :)
    endif
    
    ! surface
    ghfh (iprof,nlevels+1) = 0.

    ! medium
    Do ilayer = 2, nlevels
      z1 = gh (iprof,ilayer) ! bottom
      z2 = gh (iprof,ilayer-1) ! top
      If (z1 >= z2) then
        errorstatus = errorstatus_fatal
        Write( errMessage, '( "rttov_iniscatt_zvertinho : problem with user-defined Geopotential Height layering")' )
        Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
        if (lhook) call dr_hook('RTTOV_INISCATT_ZVERTINHO',1_jpim,zhook_handle)
        Return
      End If
      ghfh(iprof, ilayer) = (z1 + z2) / 2.
    Enddo

    ! cosmos
    ghfh (iprof,1) = 2 * gh (iprof,1) - ghfh (iprof,2)

    if (iprof .eq. 1) then
      PRINT *, "GHFH", ghfh(1, :)
    endif

  Enddo

  ! Xie: get presf 
  Do iprof = 1, nprofiles
     Do ilayer = 1, nlevels
        If (profiles (iprof) % p (ilayer) >= pressure_top) Then
            If (ilayer .ne. 1) then
              presf (iprof,ilayer) = profiles(iprof) % p(ilayer) &
              & * exp( - (ghfh (iprof, ilayer) - gh (iprof, ilayer))   / dp2dz / profiles(iprof) % t(nlevels) )
            Else
              presf (iprof,ilayer) = profiles(iprof) % p(ilayer+1) &
              & * exp( + (gh (iprof, ilayer) - ghfh (iprof, ilayer+1)) / dp2dz / profiles(iprof) % t(nlevels+1) )
            EndIf
        Else
            presf (iprof,ilayer) = pressure_top
        EndIf
    Enddo

    if (iprof .eq. 1) then
      PRINT *, "presf", presf (iprof, :) ! (nlevels)
      PRINT *, "p", profiles(iprof) % p (:) ! (nlevels+1)
    endif

  Enddo

!* Geometric variables
  Call rttov_setgeometry ( &
    & opts,            & ! in
    & .FALSE._jplm,    & ! in (inc. solar radiation)
    & .FALSE._jplm,    & ! in (plane-parallel)
    & profiles,        & ! in
    & coef=coef_rttov, & ! in
    & angles=angles)     ! out

  !* Temperature at layer boundaries (K)
  Do ichan = 1, nchannels
    iprof = chanprof(ichan) % prof
    scatt_aux % tbd (ichan,nlevels+1) = profiles (iprof) % s2m % t
    scatt_aux % tsfc (ichan)          = profiles (iprof) % skin % t
    scatt_aux % tcosmic (ichan)       = tcosmic
  Enddo

  do ilayer = 1, nlevels
    do iprof = 1, nprofiles
      tbd_prof(iprof) = profiles(iprof) % t(ilayer)
    enddo
    do ichan = 1, nchannels
      iprof = chanprof(ichan) % prof
      scatt_aux % tbd (ichan,ilayer) = tbd_prof (iprof)
    enddo
  enddo

  ! Convert temperatures to "effective" - i.e. apply band corrections
  if(coef_rttov%ff_val_bc .and. opts%rt_mw%apply_band_correction) then
    do ichan = 1, nchannels
      ichanid = chanprof(ichan) % chan
      scatt_aux%tbd(ichan,:)   = coef_rttov%ff_bco(ichanid) + coef_rttov%ff_bcs(ichanid) * scatt_aux%tbd(ichan,:)
      scatt_aux%tsfc(ichan)    = coef_rttov%ff_bco(ichanid) + coef_rttov%ff_bcs(ichanid) * scatt_aux%tsfc(ichan)    
      scatt_aux%tcosmic(ichan) = coef_rttov%ff_bco(ichanid) + coef_rttov%ff_bcs(ichanid) * scatt_aux%tcosmic(ichan)           
    enddo
  endif
  
  !* Nadir heights (km)
  Do iprof = 1, nprofiles
    Do ilayer = 1, nlevels
        scatt_aux % dz (iprof,ilayer) = ghfh (iprof,ilayer) - ghfh (iprof,ilayer+1)
     Enddo
  Enddo

  !* Get single-layer optical depths (at nadir and in hPa-1) and put onto model half levels
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof

     ! Top RTTOV level to space    
     od_rttov (1)         = -1.0_jprb * log( transmission % tau_levels (1,ichan) )     
     Do ilayer = 2, nlevels + 1 
        od_rttov (ilayer) = log( transmission % tau_levels (ilayer-1,ichan) ) &
                        & - log( transmission % tau_levels (ilayer,ichan) )
     Enddo
     ! Surface to bottom RTTOV (full pressure) level
     ! Attension: odrttov(nlevels+2) can be negative since the bottom most level can be underearth
     od_rttov (nlevels+2) = log( transmission % tau_levels (nlevels+1,ichan) ) &
                        & - log( transmission % tau_total (ichan) )

     ! Re-allocate optical depths between half pressure levels        
     od (ichan,1)         = od_rttov(1) + od_rttov(2)     
     Do ilayer = 2, nlevels - 1  
        od (ichan,ilayer) = od_rttov(ilayer+1)
     Enddo
     od (ichan,nlevels)   = od_rttov(nlevels+1) + od_rttov(nlevels+2)

     if (ichan .eq. 1) then
      PRINT *, "tau_levels", transmission % tau_levels (:, ichan)
      PRINT *, "tau_total",  transmission % tau_total (ichan)
      PRINT *, "od_rttov", od_rttov (:) ! (nlevels+2)
      PRINT *, "od", od(ichan, :) ! (nlevels)
    endif

  Enddo

  !* Change units
  Do ilayer = 1, nlevels

     !* Optical depths in km-1 and at nadir
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
     
        scatt_aux % ext (ichan,ilayer) = od (ichan,ilayer) &
          & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen 

        if (scatt_aux % ext (ichan,ilayer) < 1.0E-10_JPRB) scatt_aux % ext (ichan,ilayer) = 1.0E-10_JPRB
     Enddo
  Enddo

  !* Creation of the hydrometeor profiles
  Call rttov_hydro (        &
       & nlevels,           &! in
       & nprofiles,         &! in   
       & usercfrac,         &! in
       & presf,             &! in 
       & profiles,          &! in
       & cld_profiles,      &! in
       & coef_scatt_sets(1),&! in
       & scatt_aux)          ! inout 

  !* Cloud/rain absorption/scattering parameters
  ! WRITE(*,*) "before enter rttov_mieproc_vertinho"
  Call rttov_mieproc_vertinho (      &
       & nlevels,           &! in
       & nchannels,         &! in
       & nprofiles,         &! in
       & frequencies,       &! in
       & chanprof%prof,     &! in
       & profiles,          &! in
       & nmietables,        &! in
       & coef_scatt_sets,   &! in
       & lshapetable,       &! in
       & scatt_aux)          ! inout 
  ! WRITE(*,*) "after enter rttov_mieproc_vertinho"

  !* Scattering parameters for Eddington RT
  Call rttov_iniedd (       &
       & lradiance,         &! in
       & nlevels,           &! in
       & nchannels ,        &! in
       & nprofiles ,        &! in
       & chanprof  ,        &! in
       & angles    ,        &! in
       & coef_rttov,        &! in
       & scatt_aux)          ! inout 

  !* Surface emissivities
  zod_up_cld (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     Do ilayer = 1, nlevels
        zod_up_cld (ichan) = zod_up_cld (ichan) + scatt_aux % ext (ichan,ilayer) * scatt_aux % dz (iprof,ilayer)
     Enddo
     
     if (zod_up_cld (ichan) >= max_scatt_optical_depth) zod_up_cld (ichan) = max_scatt_optical_depth
     transmissioncld % thermal_path1 % tau_surf (0,ichan) = Exp(-1.0_JPRB * zod_up_cld (ichan) / angles (iprof) % coszen)
  Enddo
  
  Call rttov_calcemis_mw (     &
       & opts,                 &! in
       & profiles,             &! in
       & angles,               &! in
       & coef_rttov,           &! in
       & chanprof,             &! in
       & transmissioncld,      &! in
       & calcemiss,            &! in
       & scatt_aux % ems_cld,  &! inout
       & scatt_aux % ref_cld,  &! out
       & errorstatus          ) ! inout 
  
  !* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux % ems_bnd (:) = scatt_aux % ems_cld (:)
  scatt_aux % ref_bnd (:) = scatt_aux % ref_cld (:)

  !* Deallocate
  deallocate (transmissioncld % thermal_path1 % tau_surf)
  deallocate (transmissioncld % thermal_path1)

  if (lhook) call dr_hook('RTTOV_INISCATT_ZVERTINHO',1_jpim,zhook_handle)

End Subroutine rttov_iniscatt_zvertinho
