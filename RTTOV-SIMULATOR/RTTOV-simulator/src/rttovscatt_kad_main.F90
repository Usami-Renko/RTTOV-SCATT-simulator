PROGRAM rttovscatt_kad_main

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name,           &
         gas_unit_specconc

  ! rttov_types contains definitions of all RTTOV data types
  USE rttov_types, ONLY :     &
         rttov_options,       &
         rttov_options_scatt, &
         rttov_coefs,         &
         rttov_scatt_coef,    &
         rttov_profile,       &
         rttov_profile_cloud, &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  ! my mod
  USE rttovscatt_mod, ONLY : plevels, delta

  IMPLICIT NONE

#include "rttov_scatt_ad.interface"
#include "rttov_parallel_scatt_ad.interface"

#include "rttov_scatt.interface"
#include "rttov_parallel_scatt.interface"
#include "rttov_alloc_k.interface"

#include "rttov_alloc_scatt_prof.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"

#include "rttov_init_prof.interface"
#include "rttov_init_scatt_prof.interface"
#include "rttov_init_rad.interface"

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_print_opts_scatt.interface"
#include "rttov_print_profile.interface"
#include "rttov_print_cld_profile.interface"
#include "rttov_print_opts.interface"
#include "rttov_skipcommentline.interface"

! my include
#include "rttovscatt_sva.interface"
#include "rttovscatt_select_level.interface"
#include "rttovscatt_select_levels.interface"
#include "rttovscatt_gribapi.interface"
#include "rttovscatt_comput.interface"
#include "print_1d_data.interface"


  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: ioin   = 20    ! unit for input sva file
  INTEGER(KIND=jpim), PARAMETER :: ioout  = 21    ! unit for output
  INTEGER(KIND=jpim), PARAMETER :: ioouttest  = 22    ! unit for output test profile
  INTEGER(KIND=jpim), PARAMETER :: iooutMat  = 23    ! unit for output Matrix profile

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)                :: opts                          ! Options structure - leave this set to defaults
  TYPE(rttov_options_scatt)          :: opts_scatt                    ! RTTOV-SCATT options structure
  TYPE(rttov_coefs)                  :: coefs                         ! Coefficients structure
  TYPE(rttov_scatt_coef)             :: coefs_scatt                   ! RTTOV-SCATT coefficients structure
  TYPE(rttov_chanprof),      POINTER :: chanprof(:)         => NULL() ! Input channel/profile list
  INTEGER(KIND=jpim),        POINTER :: frequencies(:)      => NULL() ! Channel indexes for Mietable lookup
  LOGICAL(KIND=jplm),        POINTER :: use_chan(:,:)       => NULL() ! Flags to specify channels to simulate
  LOGICAL(KIND=jplm),        POINTER :: calcemis(:)         => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),    POINTER :: emissivity(:)       => NULL() ! Input/output surface emissivity
  TYPE(rttov_emissivity),    POINTER :: emissivity_kad(:)   => NULL() ! Input/output surface emissivity
  TYPE(rttov_profile),       POINTER :: profiles(:)         => NULL() ! Input profiles
  TYPE(rttov_profile),       POINTER :: profiles_kad(:)     => NULL() ! Input profiles
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles(:)     => NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles_kad(:) => NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles
  TYPE(rttov_radiance)               :: radiance                      ! Output radiances
  TYPE(rttov_radiance)               :: radiance_kad                     ! Output radiances

  INTEGER(KIND=jpim)                 :: errorstatus                   ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status

  ! my variables
  real(KIND=jprb), allocatable    :: levels_raw(:, :)           !(nlevels, npoints)
  real(KIND=jprb), allocatable    :: levels_selected_raw(:, :)  !(nlevels, nprof)
  real(KIND=jprb), allocatable    :: level_raw(:)                 !(npoints)
  real(KIND=jprb), allocatable    :: level_selected_raw(:)        !(nprof)

  real(KIND=jprb), allocatable    :: sva_values(:, :)

  integer(KIND=jpim)              :: nw_lat
  integer(KIND=jpim)              :: nw_lon
  integer(KIND=jpim)              :: se_lat
  integer(KIND=jpim)              :: se_lon

  real(KIND=jprb), dimension(:), allocatable    :: lats_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lons_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lats_selected_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lons_selected_raw ! (npoints)

  real(KIND=jprb), dimension(:), allocatable    :: ph
  CHARACTER(LEN=32)                             :: shortName
  real(KIND=jprb)                               :: ph_nlevel


  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: mietable_filename
  CHARACTER(LEN=256) :: sva_filename
  CHARACTER(LEN=256) :: model_filename
  CHARACTER(LEN=256) :: output_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: totalice, snowrain_units
  LOGICAL(KIND=jplm) :: use_totalice, mmr_snowrain
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  CHARACTER(LEN=11)  :: gas_unit
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch, l
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV-SCATT are as follows:
  !   1. Specify required RTTOV-SCATT options
  !   2. Read coefficients and Mietable file
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof and frequencies arrays by calling rttov_scatt_setupindex
  !   5. Read input profile(s)
  !   6. Set up surface emissivity
  !   7. Call rttov_scatt and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV-SCATT interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
  WRITE(0,*) 'enter path of Mietable file'
  READ(*,*) mietable_filename
  WRITE(0,*) 'enter path of SVA file'
  READ(*,*) sva_filename
  WRITE(0,*) 'enter path of MODEL file'
  READ(*,*) model_filename
  WRITE(0,*) 'enter path of OUTPUT file'
  READ(*,*) output_filename
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'use totalice? (0=no, 1=yes)'
  READ(*,*) totalice
  WRITE(0,*) 'snow/rain units? (0=kg/m2/s, 1=kg/kg)'
  READ(*,*) snowrain_units
  WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)
  WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads
  WRITE(0,*) 'enter nw lat'
  READ(*,*) nw_lat
  WRITE(0,*) 'enter nw lon'
  READ(*,*) nw_lon
  WRITE(0,*) 'enter se lat'
  READ(*,*) se_lat
  WRITE(0,*) 'enter se lon'
  READ(*,*) se_lon


  use_totalice = (totalice /= 0_jpim)
  mmr_snowrain = (snowrain_units /= 0_jpim)

  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV-SCATT options structure
  ! --------------------------------------------------------------------------

  ! The rttov_options structure (opts) should be left with its default values.
  ! RTTOV-SCATT only allows access to a limited number of RTTOV options: these
  ! are set in the rttov_options_scatt structure (opts_scatt).

  ! For example:
  opts_scatt % interp_mode = 1                    ! Set interpolation method
  opts_scatt % config % verbose = .False.         ! Enable printing of warnings
  opts_scatt % lradiance = .True.                 ! Use radiance weighted
  opts % config % do_checkinput = .True.          ! Disable the warning


  ! See user guide for full list of RTTOV-SCATT options


  !========== Interactive inputs == end ==============
  !===================================================


  ! --------------------------------------------------------------------------
  ! 2. Read coefficients
  ! --------------------------------------------------------------------------
  CALL rttov_read_coefs(errorstatus, coefs, opts, file_coef=coef_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! Ensure input number of channels is not higher than number stored in coefficient file
  IF (nchannels > coefs % coef % fmv_chn) THEN
    nchannels = coefs % coef % fmv_chn
  ENDIF

  ! Read the RTTOV-SCATT Mietable file
  CALL rttov_read_scattcoeffs(errorstatus, opts_scatt, coefs, coefs_scatt, file_coef=mietable_filename)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'fatal error reading RTTOV-SCATT coefficients'
    CALL rttov_exit(errorstatus)
  ENDIF

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for RTTOV direct model
  CALL rttov_alloc_k( &
        errorstatus,                &
        1_jpim,                     &  ! 1 => allocate
        nprof,                      &
        nchanprof,                  &
        nlevels,                    &
        chanprof,                   &
        opts,                       &
        profiles,                   &
        profiles_kad,               &
        coefs,                      &
        radiance = radiance,        &
        radiance_k = radiance_kad,  &
        calcemis = calcemis,        &
        emissivity = emissivity,    &
        emissivity_k = emissivity_kad,&
        init = .TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

!----------------------cld_profiles--------------------------------
  ! Allocate the RTTOV-SCATT cloud profiles structure
  ALLOCATE(cld_profiles(nprof), stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'allocation error for cld_profiles array'
    errorstatus = errorstatus_fatal
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_alloc_scatt_prof(   &
        errorstatus,             &
        nprof,                   &
        cld_profiles,            &
        nlevels,                 &
        use_totalice,            &    ! false => separate ciw and snow; true => totalice
        1_jpim,                  &    ! 1 => allocate
        init = .TRUE._jplm,      &
        mmr_snowrain = mmr_snowrain)  ! snow/rain input units: false => kg/m2/s; true => kg/kg
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for cld_profiles structure'
    CALL rttov_exit(errorstatus)
  ENDIF

!----------------------cld_profiles_kad--------------------------------
  ! Allocate the RTTOV-SCATT kad cloud profiles structure
  ALLOCATE(cld_profiles_kad(nchanprof), stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'allocation error for cld_profiles array'
    errorstatus = errorstatus_fatal
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_alloc_scatt_prof(   &
        errorstatus,             &
        nchanprof,               &
        cld_profiles_kad,        &
        nlevels,                 &
        use_totalice,            &    ! false => separate ciw and snow; true => totalice
        1_jpim,                  &    ! 1 => allocate
        init = .TRUE._jplm,      &
        mmr_snowrain = mmr_snowrain)  ! snow/rain input units: false => kg/m2/s; true => kg/kg
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for cld_profiles structure'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Populate chanprof and frequencies arrays
  ! --------------------------------------------------------------------------

  ! RTTOV-SCATT requires the frequencies array to be populated by a call to
  ! rttov_scatt_setupindex. This also populates the chanprof array. To specify
  ! only a subset of channels (i.e. those in channel_list) an array of flags is
  ! passed in (use_chan).

  ! use_chan array is dimensioned by the total number of instrument channels
  ALLOCATE(use_chan(nprof,coefs%coef%fmv_chn), &
           frequencies(nchanprof))

  ! Set use_chan to .TRUE. only for required channels
  use_chan(:,:) = .FALSE._jplm
  DO j = 1, nprof
    use_chan(j,channel_list(1:nchannels)) = .TRUE._jplm
  ENDDO

  ! Populate chanprof and frequencies arrays
  CALL rttov_scatt_setupindex ( &
        nprof,              &
        coefs%coef%fmv_chn, &
        coefs,              &
        nchanprof,          &
        chanprof,           &
        frequencies,        &
        use_chan)


  write(*, *) 'start reading profile!'
  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  profiles(:) % gas_units = 1 ! gas_units 1->kg/kg over moist air

  !===============================================
  !========== call mwri_gribapi == start =============

  !---------------[A]  levels field-----------------------

  shortName = "t"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    profiles(iprof)%t(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "q"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    profiles(iprof)%q(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "ccl"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    cld_profiles(iprof)%cc(:) = levels_selected_raw(:, iprof)/100.0
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "clwmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    cld_profiles(iprof)%clw(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "icmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    cld_profiles(iprof)%ciw(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "rwmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    cld_profiles(iprof)%rain(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "snmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof=1, nprof
    cld_profiles(iprof)%sp(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  !---------------[B] surface field--------------------

  ph_nlevel = sqrt(plevels(nlevels)*plevels(nlevels-1))

  shortName = "sp"
  call read_surface_field_latlon(model_filename, shortName, 0, level_raw, lats_raw, lons_raw)
  level_raw = level_raw / 100.0
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  call select_values(lats_raw, lats_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  call select_values(lons_raw, lons_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        IF (level_selected_raw(iprof) .gt. ph_nlevel) THEN
            profiles(iprof) % s2m % p  = level_selected_raw(iprof)
            cld_profiles(iprof) % ph(nlevels+1) = level_selected_raw(iprof)
        ELSE
            profiles(iprof) % s2m % p  = ph_nlevel + delta
            cld_profiles(iprof) % ph(nlevels+1) = ph_nlevel + delta
        ENDIF
        profiles(iprof) % latitude  =  lats_selected_raw(iprof) 
        profiles(iprof) % longitude =  lons_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)
  deallocate(lats_raw)
  deallocate(lons_raw)
  deallocate(lats_selected_raw)
  deallocate(lons_selected_raw)

  shortName = "2t"
  call read_surface_field(model_filename, shortName, 2, level_raw)
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % t  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  shortName = "q"
  call read_surface_field(model_filename, shortName, 2, level_raw)
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % q  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  shortName = "10u"
  call read_surface_field(model_filename, shortName, 10, level_raw)
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % u  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  shortName = "10v"
  call read_surface_field(model_filename, shortName, 10, level_raw)
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % v  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  shortName = "t"
  call read_surface_field(model_filename, shortName, 0, level_raw)
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        profiles(iprof) % skin % t  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  !---------------[C] p & ph--------------------

  allocate(ph(nlevels)) ! ph(1:nlevels)

  call generate_std_levels(plevels, ph)

  DO iprof = 1, nprof
      profiles(iprof) % p(:) = plevels
      cld_profiles(iprof) % ph(1:nlevels) = ph 
  ENDDO

  deallocate(ph)


  !---------------[D]  internally generated--------------

  DO iprof = 1, nprof
        profiles(iprof) % skin % salinity  = 35.0
        profiles(iprof) % skin % fastem(:) = (/3.0, 5.0, 15.0, 0.1, 0.3/)
        profiles(iprof) % skin % surftype  = 1
        profiles(iprof) % skin % watertype = 1
        profiles(iprof) % elevation = 0
  ENDDO


  !--------------[E] SVA---------------------------------

  call read_sva_fields(ioin, sva_filename, sva_values, nprof)

  DO iprof = 1, nprof
      profiles(iprof) % zenangle = sva_values(1, iprof)
      profiles(iprof) % azangle = sva_values(2, iprof)
  ENDDO

  deallocate(sva_values)

  write(*,*) "Finish reading profile!"

  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  write(*,*) "Enter rttovscatt_ad 1!"

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV-SCATT Jacobian/Adjoint model
  ! --------------------------------------------------------------------------
  
  ! The input/output K variables must be initialised to zero before every call to rttov_k:

  ! Initialise profiles_k, radiance_k and transmission_k structures to zero
  CALL rttov_init_prof(profiles_kad(:))
  CALL rttov_init_scatt_prof(cld_profiles_kad(:))
  CALL rttov_init_rad(radiance_kad)

  ! Initialise emissivity_k to zero
  emissivity_kad(:) % emis_in = 0._jprb         ! Output emissivity Jacobian is in emis_in
  emissivity_kad(:) % emis_out = 0._jprb

  ! Set input perturbation in radiance_k:
  !   If switchrad is TRUE the perturbation in bt(:) is used for "thermal" channels
  !   If switchrad is FALSE the perturbation in total(:) is used
  !   For solar-only channels the perturbation in total(:) is *always* used
  !   For rttovscatt channels the perturbation in bt(:) is *always* used
  !   It is harmless to specify inputs in both bt/total in any case: RTTOV will use the
  !     appropriate input for each channel
  radiance_kad % total(:) = 1._jprb
  radiance_kad % bt(:) = 1._jprb

  write(*,*) "Enter rttovscatt_ad 2!"

  IF (nthreads <= 1) THEN
    CALL rttov_scatt_ad ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for Mietable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coefs,               &! in    coefficients structure
        coefs_scatt,         &! in    Mietable structure
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        profiles_kad,        &! inout
        cld_profiles_kad,    &! inout
        emissivity_kad,      &! inout 
        radiance,            &! inout computed radiances
        radiance_kad)         ! inout
  ELSE
    CALL rttov_parallel_scatt_ad ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for Mietable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coefs,               &! in    coefficients structure
        coefs_scatt,         &! in    Mietable structure
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        profiles_kad,        &! inout
        cld_profiles_kad,    &! inout
        emissivity_kad,      &! inout 
        radiance,            &! inout computed radiances
        radiance_kad,        &! inout 
        nthreads = nthreads)  ! in    number of threads to use
  ENDIF
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_scatt error'
    CALL rttov_exit(errorstatus)
  ENDIF

  !---------------------------------------------------------------------------
  !  *. output the profile for a test
  !---------------------------------------------------------------------------

  OPEN(ioouttest, file='output_test_profile.dat', status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

!---------------------------1----------------------------------------------------

  WRITE(ioouttest,*)' -----------------'
  WRITE(ioouttest,*)' Instrument ', inst_name(coefs % coef % id_inst)
  WRITE(ioouttest,*)' -----------------'
  WRITE(ioouttest,*)' '
  CALL rttov_print_opts(opts, lu=ioouttest)

  iprof = 1

  joff = (iprof-1_jpim) * nchannels

  nprint = 1 + INT((nchannels-1)/10)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)' Profile ', iprof

  CALL rttov_print_profile(profiles(iprof), lu=ioouttest)

  WRITE(ioouttest,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
  WRITE(ioouttest,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
  WRITE(ioouttest,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
  WRITE(ioouttest,222) (radiance % total(j), j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED OVERCAST RADIANCES:'
  WRITE(ioouttest,222) (radiance % cloudy(j), j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED SURFACE EMISSIVITIES:'
  WRITE(ioouttest,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)

  WRITE(ioouttest,*)
  WRITE(ioouttest,*)

!------------------------------------109----------------------------------------

  CALL rttov_print_opts(opts, lu=ioouttest)

  iprof = 109

  joff = (iprof-1_jpim) * nchannels

  nprint = 1 + INT((nchannels-1)/10)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)' Profile ', iprof

  CALL rttov_print_profile(profiles(iprof), lu=ioouttest)

  WRITE(ioouttest,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
  WRITE(ioouttest,111) (chanprof(j) % chan, j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED BRIGHTNESS TEMPERATURES (K):'
  WRITE(ioouttest,222) (radiance % bt(j), j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED RADIANCES (mW/m2/sr/cm-1):'
  WRITE(ioouttest,222) (radiance % total(j), j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED OVERCAST RADIANCES:'
  WRITE(ioouttest,222) (radiance % cloudy(j), j = 1+joff, nchannels+joff)
  WRITE(ioouttest,*)' '
  WRITE(ioouttest,*)'CALCULATED SURFACE EMISSIVITIES:'
  WRITE(ioouttest,444) (emissivity(j) % emis_out, j = 1+joff, nchannels+joff)

  WRITE(ioouttest,*)
  WRITE(ioouttest,*)


  ! Close output file
  CLOSE(ioouttest, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF
  
  !=====================================================
  !============== Output results == start ==============

  ! Open output file where results are written
  OPEN(ioout, file=TRIM(output_filename), status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  ! WRITE(ioout,*)' -----------------'
  ! WRITE(ioout,*)' Instrument ', inst_name(coefs % coef % id_inst)
  ! WRITE(ioout,*)' -----------------'
  ! WRITE(ioout,*)' '
  ! CALL rttov_print_opts_scatt(opts_scatt, lu=ioout)

  !-----------------------reference---------------------------------

  WRITE(ioout,777)'CHANNELS PROCESSED FOR SAT ', platform_name(coefs % coef % id_platform), coefs % coef % id_sat
  WRITE(ioout,111) (chanprof(j) % chan, j = 1, nchannels)
  WRITE(ioout, *) ' '

  IF (profiles(1) % gas_units == gas_unit_specconc) THEN
    gas_unit = ' (K/(kg/kg)) '
  ELSE
  gas_unit = ' (K/ppmv) '
  ENDIF

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    WRITE(ioout,*)' Profile ', iprof
    CALL rttov_print_cld_profile(cld_profiles(iprof), lu=ioout)    

    WRITE(ioout,444) (radiance % bt(j), j = 1+joff, nchannels+joff)


    DO j = 1, nchannels

      WRITE(ioout,*)'Channel ', chanprof(j+joff) % chan
      WRITE(ioout,'(/a5,a9,a18,1x,a19,6x,a12,6x,a12,6x,a12,6x,a12,6x,a12)',advance='no') &
      "Level", "Pressure", "JAC(Temp) (K/K)", "JAC(WV) "//gas_unit, &
      "rain", "sp", "clw", "ciw", "cc"
      WRITE(ioout,'(a)',advance='yes')

      DO l = 1, profiles_kad(j+joff) % nlevels
        WRITE(ioout,'(i4,1x,f9.4,1x,e16.6,1x,e17.6,1x,e17.6,1x,e17.6,1x,e17.6,1x,e17.6,1x,e17.6)',advance='no') &
          l, profiles(iprof)%p(l), profiles_kad(j+joff)%t(l), profiles_kad(j+joff)%q(l), &
          cld_profiles_kad(j+joff)%rain(l), cld_profiles_kad(j+joff)%sp(l),&
          cld_profiles_kad(j+joff)%clw(l), cld_profiles_kad(j+joff)%ciw(l),&
          cld_profiles_kad(j+joff)%cc(l)
        WRITE(ioout,'(a)',advance='yes')
      ENDDO
      WRITE(ioout,*)' '

  ENDDO

  ENDDO


  ! Close output file
  CLOSE(ioout, iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error closing the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  !============== Output results == end ==============
  !=====================================================


  ! --------------------------------------------------------------------------
  ! 8. Deallocate all RTTOV arrays and structures
  ! --------------------------------------------------------------------------
  DEALLOCATE(channel_list, use_chan, frequencies, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  CALL rttov_alloc_scatt_prof(   &
        errorstatus,             &
        nprof,                   &
        cld_profiles,            &
        nlevels,                 &
        use_totalice,            &  ! must match value used for allocation
        0_jpim)                     ! 0 => deallocate
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for cld_profiles structure'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_alloc_scatt_prof(   &
        errorstatus,             &
        nchanprof,               &
        cld_profiles_kad,        &
        nlevels,                 &
        use_totalice,            &  ! must match value used for allocation
        0_jpim)                     ! 0 => deallocate
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for cld_profiles_kad structure'
    CALL rttov_exit(errorstatus)
  ENDIF

  DEALLOCATE(cld_profiles, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'dellocation error for cld_profiles array'
  ENDIF

  DEALLOCATE(cld_profiles_kad, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'dellocation error for cld_profiles_kad array'
  ENDIF

  ! Deallocate structures for RTTOV Jacobian model
  CALL rttov_alloc_k( &
        errorstatus,                &
        0_jpim,                     &  ! 0 => deallocate
        nprof,                      &
        nchanprof,                  &
        nlevels,                    &
        chanprof,                   &
        opts,                       &
        profiles,                   &
        profiles_kad,               &
        coefs,                      &
        radiance = radiance,        &
        radiance_k = radiance_kad,    &
        calcemis = calcemis,        &
        emissivity = emissivity,    &
        emissivity_k = emissivity_kad,&
        init = .TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_scattcoeffs(coefs_scatt)

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A9,I3)

END PROGRAM rttovscatt_kad_main
