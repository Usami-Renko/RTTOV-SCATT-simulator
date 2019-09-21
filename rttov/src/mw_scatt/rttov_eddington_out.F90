!
Subroutine rttov_eddington_out ( &
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & chanprof,          &! in
     & angles,            &! in
     & scatt_aux,         &! in
     & rad_cld,           &! out 
     & sfc_terms)          ! out, optional, ! Downward radiance source terms, Upward radiance source terms, Total transmittancs

  ! Description:
  ! to compute Eddington approximation to RT
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
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and 
  !     precipitation. Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer 
  !     modeling in clouds and precipitation. Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       09/2002   Initial version     (P. Bauer, E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV-9 version  / includes linear in tau
  !                      approximation for the clear sky source function (A. Geer)
  !  1.6       07/2008   Speed up and make consistent minimum SSA (A. Geer)
  !  1.7       07/2008   Clear sky speed-ups (A. Geer)
  !  1.8       06/2009   Bug fixed in index (T. Wilhelmsson)
  !  1.9       03/2010   Use mclayer rather than min_ssa (A. Geer)
  !            02/2013   Return radiance source terms (OPTIONAL): downward&upward radiance, total transmittance (F. Baordo)
  !            11/2017   R/T now done with radiances, not Tb (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:
  
  Use rttov_types, Only :      &
       & rttov_chanprof       ,&
       & rttov_geometry       ,&
       & rttov_profile_scatt_aux ,&
       & eddington_sfc_type

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  Use rttov_const, Only: ccthres

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of channels*profiles=radiances
  type(rttov_chanprof), intent (in) :: chanprof(nchannels)  ! Channel and profile indices
  
  Type (rttov_geometry),            Intent (in)  :: angles    (nprofiles) ! Zenith angles
  Type (rttov_profile_scatt_aux),   Intent (in)  :: scatt_aux             ! Auxiliary profile variables for RTTOV_SCATT
  Real (Kind=jprb),                 Intent (out) :: rad_cld   (nchannels) ! Radiances

  ! Downward radiance source terms, Upward radiance source terms, Total transmittances
  Type (eddington_sfc_type), optional, Intent (inout) :: sfc_terms       
      
!INTF_END

!* Local variables
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dp   ! D+ for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: dm   ! D- for boundary conditions
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_up ! Upward radiance source terms
  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_do ! Downward radiance source terms

  Real (Kind=jprb), Dimension (nchannels) :: irad_do, ftop ! Downward radiances
  Real (Kind=jprb), Dimension (nchannels) :: irad_up       ! Upward radiances
  Real (Kind=jprb), Dimension (nchannels) :: irad_sfc      ! Inward radiances at surface
  Real (Kind=jprb), Dimension (nchannels) :: tau_t         ! Total transmittancs

  Real (Kind=jprb), Dimension (nchannels,nlevels) :: j_part  ! Part of the linear in tau

  Integer (Kind=jpim)  :: ilayer, jlayer, iprof, ichan
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

#include "rttov_boundaryconditions.interface"
#include "rttov_integratesource.interface"

  IF (LHOOK) CALL DR_HOOK('RTTOV_EDDINGTON_OUT',0_jpim,ZHOOK_HANDLE)        

  j_up (:,:) = 0.0_JPRB
  j_do (:,:) = 0.0_JPRB

!* Channels * Profiles      
  do ichan = 1, nchannels
    iprof   = chanprof(ichan)%prof
    
    !* Top
    irad_sfc (ichan) = scatt_aux % ems_cld (ichan) * scatt_aux % bsfc (ichan)
    
    if (scatt_aux % cfrac (iprof) > ccthres) then

      !* Clear-sky source terms (but only used outside Eddington)
      !* Linear in tau for consistency with RTTOV9
      do ilayer = 1, nlevels
        if (ilayer < scatt_aux % mclayer(ichan)) then
          if ( scatt_aux % delta (ichan,ilayer) > 1E-8 ) then

            j_part (ichan,ilayer) = &
             &   (scatt_aux % b0 (ichan,ilayer) - scatt_aux % bn (ichan,ilayer)) &
             & * (1.0_JPRB  - scatt_aux % tau (ichan,ilayer)) &
             & / scatt_aux % delta (ichan,ilayer)

            j_up (ichan,ilayer)   = scatt_aux % bn (ichan,ilayer) &
             & - scatt_aux % b0 (ichan,ilayer) * scatt_aux % tau (ichan,ilayer) &
             & + j_part(ichan,ilayer)

            j_do (ichan,ilayer)   = scatt_aux % b0 (ichan,ilayer) &
             & - scatt_aux % bn (ichan,ilayer) * scatt_aux % tau (ichan,ilayer) &
             & - j_part(ichan,ilayer)

          else

            ! Linear in tau formula is computationally unstable for small delta - avoid
            j_up (ichan,ilayer) = scatt_aux % b0 (ichan,ilayer) &
             &                  * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 
            j_do (ichan,ilayer) = scatt_aux % b0 (ichan,ilayer) &
             &                  * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 

          endif
        endif
      enddo

      !* Downward radiance at cloud top
      irad_do (ichan) = scatt_aux % btop (ichan)

      Do ilayer = 1, scatt_aux % mclayer (ichan) - 1
        irad_do (ichan) = irad_do (ichan) * scatt_aux % tau (ichan,ilayer) &
          &             + j_do (ichan,ilayer)
      End do
 
      ftop (ichan) = irad_do (ichan)

    endif
  End do

!* Get D+, D- from boundary conditions
     Call rttov_boundaryconditions (&
&     nlevels,       &! in
&     nchannels,     &! in
&     chanprof%prof, &! in
&     scatt_aux,     &! in
&     ftop,          &! in
&     dp,            &! out
&     dm)             ! out 

!* Integrate radiance source terms
     Call rttov_integratesource (&
&     nlevels,       &! in
&     nchannels,     &! in
&     nprofiles,     &! in
&     chanprof%prof, &! in
&     angles,        &! in
&     scatt_aux,     &! in
&     dp,            &! in
&     dm,            &! in
&     j_do,          &! inout
&     j_up)           ! inout 

!* Integrate downward radiances/transmittance
  irad_do (:) = scatt_aux % btop (:)
  irad_up (:) = irad_sfc   (:)
  tau_t   (:) = 1.0_JPRB
  
  Do ilayer = 1, nlevels
     jlayer = nlevels + 1 - ilayer
  
     irad_do (:) = irad_do (:) * scatt_aux % tau (:,ilayer) + j_do (:,ilayer)
     irad_up (:) = irad_up (:) * scatt_aux % tau (:,jlayer) + j_up (:,jlayer)
     
     tau_t   (:) = tau_t   (:) * scatt_aux % tau (:,jlayer)
  Enddo

  rad_cld (:) = irad_up (:) + scatt_aux % ref_cld (:) * irad_do (:) * tau_t (:)
    
  if(present(sfc_terms)) then
    sfc_terms % down = irad_do (:)
    sfc_terms % up   = irad_up (:)
    sfc_terms % tau  = tau_t   (:)
  endif
  

  IF (LHOOK) CALL DR_HOOK('RTTOV_EDDINGTON_OUT',1_jpim,ZHOOK_HANDLE)

End Subroutine rttov_eddington_out
