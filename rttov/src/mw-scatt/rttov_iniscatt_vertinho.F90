!
Subroutine rttov_iniscatt_vertinho (&
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
      & scatt_aux)          ! inout 

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

  Use parkind1, Only : jpim, jplm
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

!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, ichan, ichanid
  Real    (Kind=jprb) :: p1, p2, pm, dp2dz

  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf  ! Pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh ! Half-level pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop  ! Thickness of top half of level [km]
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot  ! Thickness of bottom half of level [km]
  Real (Kind=jprb), Dimension (2:nlevels)           :: dzr    ! Thickness of RTTOV level [km]
  Real (Kind=jprb), Dimension (nchannels,nlevels)   :: od     ! Single layer optical depths on our levels
  Real (Kind=jprb), Dimension (nlevels+1)           :: od_rttov ! Single layer optical depths on RTTOV levels
  Real (Kind=jprb), Dimension (nchannels)           :: zod_up_cld   ! Optical depth from top of the atmosphere 
  Real (Kind=jprb), Dimension (nprofiles)           :: tbd_prof ! Boundary temperature, by profile
  
  Type (rttov_transmission_aux) :: transmissioncld  ! Clear+cloud transmittances with cloud

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_iniscatt '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
#include "rttov_mieproc_vertinho.interface"
#include "rttov_iniedd.interface"
#include "rttov_hydro.interface"
#include "rttov_calcemis_mw.interface"
#include "rttov_setgeometry.interface"
#include "rttov_errorreport.interface"

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INISCATT_VERTINHO',0_jpim,zhook_handle)

  ! write(*,*) "[in]: rttov_iniscatt_vertinho"

  errorstatus = errorstatus_success

  ! NB this only truly needed if using FASTEM 3 
  allocate (transmissioncld % thermal_path1)
  allocate (transmissioncld % thermal_path1 % tau_surf (0:0,nchannels))

  dp2dz = - rgc / gravity / mair

  scatt_aux % ext (:,:) = 0.0_JPRB
  scatt_aux % ssa (:,:) = 0.0_JPRB
  scatt_aux % asm (:,:) = 0.0_JPRB

!* Security on user-defined pressures
  Do iprof = 1, nprofiles
     Do ilayer = 1, nlevels
        If (profiles (iprof) % p (ilayer) >= pressure_top) Then
            presf (iprof,ilayer) = profiles (iprof) % p  (ilayer)
        else
            presf (iprof,ilayer) = pressure_top
        endif
    Enddo
    Do ilayer = 1, nlevels + 1
       If (cld_profiles (iprof) % ph (ilayer) >= pressure_top) Then
           presfh (iprof,ilayer) = cld_profiles (iprof) % ph (ilayer)
       else
           presfh (iprof,ilayer) = pressure_top
       endif
    Enddo
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
    scatt_aux % tbd (ichan,1)         = profiles (iprof) % t (1)
    scatt_aux % tsfc (ichan)          = profiles (iprof) % skin % t
    scatt_aux % tcosmic (ichan)       = tcosmic
  Enddo

  do ilayer = 1, nlevels-1
    do iprof = 1, nprofiles
      p1 = presf  (iprof,ilayer+1)
      p2 = presf  (iprof,ilayer  )
      pm = presfh (iprof,ilayer+1)
      tbd_prof (iprof) = profiles (iprof) % t (ilayer+1) + (profiles (iprof) % t (ilayer) & 
                     & - profiles (iprof) % t (ilayer+1)) / log(p2/p1) * log(pm/p1) 
    enddo
    do ichan = 1, nchannels
      iprof = chanprof(ichan) % prof
      scatt_aux % tbd (ichan,ilayer+1) = tbd_prof (iprof)
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
  Do ilayer = nlevels, 1, -1
     Do iprof = 1, nprofiles
        p1 = presfh (iprof,ilayer+1)
        p2 = presfh (iprof,ilayer  )
        pm = presf  (iprof,ilayer  )

        If (p1 <= p2) then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "iniscatt : problem with user-defined pressure layering")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
           if (lhook) call dr_hook('RTTOV_INISCATT',1_jpim,zhook_handle)
           Return
        End If

        scatt_aux % dz (iprof,ilayer) = dp2dz * Log(p2/p1) * profiles (iprof) % t (ilayer)

        dzbot (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles (iprof) % t (ilayer)
        dztop (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles (iprof) % t (ilayer)

     Enddo
  Enddo

  !* Get single-layer optical depths (at nadir and in hPa-1) and put onto model half levels
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof

     ! Top RTTOV level to space    
     od_rttov (1)         = -1.0_jprb * log( transmission % tau_levels (1,ichan) )     
     Do ilayer = 2, nlevels 
        od_rttov (ilayer) = log( transmission % tau_levels (ilayer-1,ichan) ) &
                        & - log( transmission % tau_levels (ilayer,ichan) )
     Enddo
     ! Surface to bottom RTTOV (full pressure) level
     od_rttov (nlevels+1) = log( transmission % tau_levels (nlevels,ichan) ) &
                        & - log( transmission % tau_total (ichan) )

     ! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzr(ilayer) = dzbot(iprof,ilayer-1) + dztop(iprof,ilayer)
     Enddo

     ! Re-allocate optical depths between half pressure levels        
     od (ichan,1)         = od_rttov(1) &
                        & + od_rttov(2) * dzbot(iprof,1) / dzr(2)     
     Do ilayer = 2, nlevels - 1  
        od (ichan,ilayer) = od_rttov(ilayer)   * dztop(iprof,ilayer) / dzr(ilayer) &
                        & + od_rttov(ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) 
     Enddo
     od (ichan,nlevels)   = od_rttov(nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + od_rttov(nlevels+1) 

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

  if (lhook) call dr_hook('RTTOV_INISCATT_VERTINHO',1_jpim,zhook_handle)

End Subroutine rttov_iniscatt_vertinho
