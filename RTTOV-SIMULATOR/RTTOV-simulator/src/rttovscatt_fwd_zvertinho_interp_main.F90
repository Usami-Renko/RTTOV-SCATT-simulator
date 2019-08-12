PROGRAM rttovscatt_fwd_zvertinho_interp_main

  ! rttov_const contains useful RTTOV constants
  USE rttov_const, ONLY :     &
         errorstatus_success, &
         errorstatus_fatal,   &
         platform_name,       &
         inst_name

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
  
  USE rttovscatt_mod, ONLY : model_res, base_nw_lon, base_nw_lat, Ni, Nj, &
  & rgc, mair, gravity

  IMPLICIT NONE

#include "rttov_scatt.interface"
#include "rttov_scatt_zvertinho.interface"
#include "rttov_parallel_scatt.interface"
#include "rttov_alloc_scatt_prof.interface"
#include "rttov_read_scattcoeffs.interface"
#include "rttov_dealloc_scattcoeffs.interface"
#include "rttov_scatt_setupindex.interface"

#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_direct.interface"
#include "rttov_print_opts_scatt.interface"
#include "rttov_print_profile.interface"
#include "rttov_print_cld_profile.interface"
#include "rttov_skipcommentline.interface"

! my include
#include "rttovscatt_sva.interface"
#include "rttovscatt_select_level.interface"
#include "rttovscatt_select_levels.interface"
#include "rttovscatt_readbinary.interface"
#include "rttovscatt_comput.interface"
#include "rttov_print_opts.interface"


  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: ioin   = 20    ! unit for input sva file
  INTEGER(KIND=jpim), PARAMETER :: iomd   = 89    ! unit for binary modelvar
  INTEGER(KIND=jpim), PARAMETER :: ioout  = 21    ! unit for output
  INTEGER(KIND=jpim), PARAMETER :: ioouttest  = 22    ! unit for output test profile

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)                :: opts                     ! Options structure - leave this set to defaults
  TYPE(rttov_options_scatt)          :: opts_scatt               ! RTTOV-SCATT options structure 
  TYPE(rttov_coefs)                  :: coefs                    ! Coefficients structure
  TYPE(rttov_chanprof),      POINTER :: chanprof(:)    => NULL() ! Input channel/profile list
  INTEGER(KIND=jpim),        POINTER :: frequencies(:) => NULL() ! Channel indexes for Mietable lookup
  LOGICAL(KIND=jplm),        POINTER :: use_chan(:,:)  => NULL() ! Flags to specify channels to simulate
  LOGICAL(KIND=jplm),        POINTER :: calcemis(:)    => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),    POINTER :: emissivity(:)  => NULL() ! Input/output surface emissivity
  TYPE(rttov_profile),       POINTER :: profiles(:)    => NULL() ! Input profiles
  TYPE(rttov_profile_cloud), POINTER :: cld_profiles(:)=> NULL() ! Input RTTOV-SCATT cloud/hydrometeor profiles
  TYPE(rttov_radiance)               :: radiance                 ! Output radiances

  INTEGER(KIND=jpim)                 :: errorstatus              ! Return error status of RTTOV subroutine calls

  INTEGER(KIND=jpim) :: alloc_status

  ! my variables
  real(KIND=jprb), allocatable    :: levels_raw(:, :)           !(nlevels, npoints)
  real(KIND=jprb), allocatable    :: levels_raw2(:, :)          !(nlevels, npoints)
  real(KIND=jprb), allocatable    :: levels_selected_raw(:, :)  !(nlevels, nprof)
  real(KIND=jprb), allocatable    :: levels_selected_raw2(:, :) !(nlevels, nprof) 
  real(KIND=jprb), allocatable    :: level_raw(:)                 !(npoints)
  real(KIND=jprb), allocatable    :: level_selected_raw(:)        !(nprof)

  ! interp module modified
  real(KIND=jprb), allocatable    :: sva_values(:, :)
  integer(KIND=jpim), allocatable :: select_table(:, :)

  real(KIND=jprb), dimension(:, :), allocatable :: lats
  real(KIND=jprb), dimension(:, :), allocatable :: lons

  real(KIND=jprb), dimension(:), allocatable    :: lats_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lons_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lats_selected_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lons_selected_raw ! (npoints)

  real(KIND=jprb), dimension(:), allocatable    :: ph !(nlevels) half-layer pressure
  
  ! zvertinho variables
  integer(KIND=jpim)                            :: b_rects ! to call readbinary
  character(LEN=32)                             :: varname ! varname for test
  real(KIND=jprb), dimension(:, :), allocatable :: gh      ! geopotential height
  integer(KIND=jpim)                            :: ilon, ilat
  real(KIND=jprb), dimension(:), allocatable    :: half_container
  real(KIND=jprb)                               :: dp2dz

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: mietable_dir
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
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  ! variables for input vertical inhomogeneity
  integer(KIND=jpim) :: vertinho_mode
  integer(KIND=jpim) :: nmietables
  CHARACTER(LEN=256), dimension(:), allocatable :: mietable_filenames
  integer(KIND=jpim) :: nshapelayers
  integer(KIND=jpim), dimension(:), allocatable :: lshape
  integer(KIND=jpim), dimension(:), allocatable :: lshape_bot
  integer(KIND=jpim), dimension(:,:), allocatable :: lshapetable

  !verticle inhomogeneity loop variables
  integer(KIND=jpim) :: ishapelayer
  integer(KIND=jpim) :: imietable

  ! Mietable structure
  TYPE(rttov_scatt_coef), POINTER  :: coefs_scatt_sets(:) => NULL()              ! RTTOV-SCATT coefficients structure

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

  ! WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
  ! WRITE(0,*) 'enter path of tempSVA file'
  READ(*,*) sva_filename
  ! WRITE(0,*) 'enter path of MODEL file'
  READ(*,*) model_filename
  ! WRITE(0,*) 'enter path of OUTPUT file'
  READ(*,*) output_filename
  ! WRITE(0,*) 'enter number of valid_profiles'
  READ(*,*) nprof
  ! WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels ! Not include the underearth virtual half level
  ! WRITE(0,*) 'use totalice? (0=no, 1=yes)'
  READ(*,*) totalice
  ! WRITE(0,*) 'snow/rain units? (0=kg/m2/s, 1=kg/kg)'
  READ(*,*) snowrain_units
  ! WRITE(0,*) 'enter number of channels to simulate per profile'
  READ(*,*) nchannels
  ALLOCATE(channel_list(nchannels))
  ! WRITE(0,*) 'enter space-separated channel list'
  READ(*,*,iostat=ios) channel_list(:)
  ! WRITE(0,*) 'enter number of threads to use'
  READ(*,*) nthreads

  !input variables of vertical inhomogeneity
  ! WRITE(0,*) 'enter vertical inhomogeneity mode'
  READ(*,*) vertinho_mode
  ! WRITE(0,*) 'enter number of mietable files'
  READ(*,*) nmietables
  ! WRITE(0,*) 'enter directory of mietable'
  READ(*,*) mietable_dir

  ALLOCATE(mietable_filenames(nmietables))
  ! WRITE(0,*) 'enter path of Mietable filenames'
  READ(*,*) mietable_filenames
  DO imietable = 1, nmietables
    ! WRITE(0,*) trim(mietable_filenames(imietable))
  ENDDO

  IF(vertinho_mode .eq. 1) THEN
    ! WRITE(0,*) 'enter number of shape layers'
    READ(*,*) nshapelayers

    ALLOCATE(lshape(nshapelayers))
    ALLOCATE(lshape_bot(nshapelayers))
    
    ! WRITE(0,*) 'enter shape index of shape layers'
    READ(*,*) lshape

    ! WRITE(0,*) 'enter shape bottom level index of shape layers'
    READ(*,*) lshape_bot
  ENDIF


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
  opts_scatt % config % apply_reg_limits = .True. ! restrict to the regression limit
  opts_scatt % config % do_checkinput = .True.    ! enable the warning



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

  ! Read the RTTOV-SCATT Mietable file (nmietables)
  ALLOCATE(coefs_scatt_sets(nmietables))
  DO imietable = 1, nmietables
    mietable_filename = trim(mietable_dir)//trim(mietable_filenames(imietable))
    ! WRITE(0,*) mietable_filename

    CALL rttov_read_scattcoeffs(errorstatus, opts_scatt, coefs, & 
      & coefs_scatt_sets(imietable), file_coef=mietable_filename)
    IF (errorstatus /= errorstatus_success) THEN
      WRITE(*,*) 'fatal error reading RTTOV-SCATT coefficients'
      CALL rttov_exit(errorstatus)
    ENDIF

  ENDDO

  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for RTTOV direct model
  CALL rttov_alloc_direct( &
        errorstatus,             &
        1_jpim,                  &  ! 1 => allocate
        nprof,                   &
        nchanprof,               &
        nlevels+1,               &  ! the underearth half level
        chanprof,                &
        opts,                    &
        profiles,                &
        coefs,                   &
        radiance = radiance,     &
        calcemis = calcemis,     &
        emissivity = emissivity, &
        init = .TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

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
  !========== call gribapi == start =============

  !----------------[A0] SVA---------------------------------

  call read_unstructured_sva_fields(ioin, sva_filename, sva_values, select_table, nprof)

  DO iprof = 1, nprof
      profiles(iprof) % zenangle = sva_values(1, iprof)
      profiles(iprof) % azangle = sva_values(2, iprof)
  ENDDO

  deallocate(sva_values)

  !---------------[A1] Longitude & Latitude---------------

  allocate(lats(Ni, Nj))
  allocate(lons(Ni, Nj))
  allocate(lats_raw(nprof))
  allocate(lons_raw(nprof))

  DO ilon = 1, Ni
    DO ilat = 1, Nj
      lats(ilon, ilat) = base_nw_lat - model_res * (ilat - 1)
      lons(ilon, ilat) = base_nw_lon + model_res * (ilon - 1)
    ENDDO
  ENDDO

  lats_raw = reshape(lats, (/Ni*Nj/))
  lons_raw = reshape(lons, (/Ni*Nj/))

  call select_unstructured_values(lats_raw, lats_selected_raw, select_table, nprof)
  call select_unstructured_values(lons_raw, lons_selected_raw, select_table, nprof)

  DO iprof = 1, nprof
    profiles(iprof) % latitude = lats_selected_raw(iprof)
    profiles(iprof) % longitude = lons_selected_raw(iprof)
  ENDDO
  
  deallocate(lats)
  deallocate(lons)
  deallocate(lats_raw)
  deallocate(lons_raw)
  deallocate(lats_selected_raw)
  deallocate(lons_selected_raw)

  !---------------[A2] original full levels--------------------

  ! allocate(ph(nlevels))  ! ph(1:nlevels)

  ! get pi(nlevels+1, nprof) on half level
  varname = "pi"
  b_rects = 50
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels+1)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)

  ! get th(nlevels, nprof) on full level
  allocate(half_container(nlevels+1))
  varname = "th"
  b_rects = 101
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw2, nlevels)
  call select_levels_unstructured_values(levels_raw2, levels_selected_raw2, select_table, nprof)

  DO iprof = 1, nprof
    
    !---------- populate profiles(:)%t(1:levels+1) = pi * th on the original half level
    ! do interpolating and extrapolating
    half_container = 0.
    half_container(2:nlevels+1) = levels_selected_raw2(:, iprof) 
    half_container(1:nlevels)   = half_container(1:nlevels) + levels_selected_raw2(:, iprof)
    half_container              = half_container / 2.
    half_container(1)           = 2 * levels_selected_raw2(2, iprof) - levels_selected_raw2(3, iprof)
    half_container(nlevels+1)   = 2 * levels_selected_raw2(nlevels, iprof) - levels_selected_raw2(nlevels-1, iprof)
    ! T = pi * theta
    profiles(iprof) % t(:) = levels_selected_raw(:, iprof) * half_container

    !---------- populate profiles(:)%p(1:nlevels)  pi --> p on the original half level
    profiles(iprof) % p(:) = 1000 * ((levels_selected_raw(:, iprof)) ** 3.5) 

    !---------- populate profiles(:)%ph(1:nlevels)
    ! seems no need for ph any more
    ! CALL generate_std_levels(profiles(iprof)%p, ph)
    ! cld_profiles(iprof)%ph(1:nlevels) = ph(:) 
    !-----------other variables
    ! profiles(iprof) % s2m % p  = 1000 * log(levels_selected_raw(nlevels, iprof))/log(0.2857) 
    ! cld_profiles(iprof) % ph(nlevels+1) = profiles(iprof) % s2m % p

    ! test
    IF (iprof .eq. 1) THEN
      write(*,*) 'p[1]=', profiles(iprof)%p(1), 'p[50]=', profiles(iprof)%p(50), 'p[51]=', profiles(iprof)%p(51)
      write(*,*) 't[1]=', profiles(iprof)%t(1), 't[50]=', profiles(iprof)%t(50), 't[51]=', profiles(iprof)%t(51)
    ENDIF
  ENDDO

  ! deallocate(ph)
  deallocate(half_container)
  deallocate(levels_raw)
  deallocate(levels_raw2)
  deallocate(levels_selected_raw)
  deallocate(levels_selected_raw2)
  
  ! get q
  b_rects = 299 ! Qv
  varname = "Qv"
  allocate(half_container(nlevels+1))
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)
  DO iprof=1, nprof
    ! do interpolating and extrapolating
    half_container = 0.
    half_container(2:nlevels+1) = levels_selected_raw(:, iprof) 
    half_container(1:nlevels)   = half_container(1:nlevels) + levels_selected_raw(:, iprof)
    half_container              = half_container / 2.
    half_container(1)           = 2 * levels_selected_raw(2, iprof) - levels_selected_raw(3, iprof)
    half_container(nlevels+1)   = 2 * levels_selected_raw(nlevels, iprof) - levels_selected_raw(nlevels-1, iprof)

    ! avoid negative q
    If (half_container(1) .le. 0.) THEN
        half_container(1) = 0. 
    Endif

    If (half_container(nlevels+1) .le. 0.) THEN
        half_container(nlevels+1) = 0.
    Endif

    profiles(iprof)%q(:) = half_container

    IF (iprof .eq. 1) THEN
      write(*,*) 'q[1]=', profiles(iprof)%q(1), 'q[50]=', profiles(iprof)%q(50), 'q[51]=', profiles(iprof)%q(51)
    ENDIF
  ENDDO
  deallocate(half_container)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  !--------------[A3] gh ------------------------------
  allocate(gh(nprof, nlevels))

  dp2dz = - rgc / gravity / mair

  ! get gh(nlevels, nprof)
  varname = "gh"
  b_rects = 0
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)

  ! get zs(nprof)  
  varname = "zs"
  b_rects = 849 !zs
  call read_surface_field(iomd, model_filename, varname, b_rects, level_raw)
  call select_unstructured_values(level_raw, level_selected_raw, select_table, nprof)

  DO iprof = 1, nprof
    gh(iprof, :) = (levels_selected_raw(:, iprof) - level_selected_raw(iprof)) / 1000. ! convert from [m] --> [km]
  ENDDO

  ! hydrostatic inference of 2mp
  DO iprof = 1, nprof 
    profiles(iprof) % s2m % p  = profiles(iprof) % p(nlevels) &
    & * exp( - gh(iprof, nlevels-1) / 2. / dp2dz / profiles(iprof) % t(nlevels) )

    IF (iprof .eq. 1) THEN
      write(*,*) 'gh[1]=', gh(iprof, 1), 'gh[50]=', gh(iprof, 50)
      write(*,*) '2mp', profiles(iprof) % s2m % p
    ENDIF

  ENDDO

  ! fix this later
  gh = gh + 0.001 ! add 1m to avoid the sealevel dzbot = 0

  deallocate(level_raw)
  deallocate(level_selected_raw)
  deallocate(levels_raw)
  deallocate(levels_selected_raw)

  !---------------[A4] populate the lshapetable---------
  ! vertical inhomogeniety mode 1 : shapelayers

  IF(vertinho_mode .eq. 1) THEN

    allocate(lshapetable(nprof, nlevels))

    DO iprof = 1, nprof

      !WRITE(*, *) "[Main]: vertinho 1"

      IF(nshapelayers .eq. 1) THEN
        lshapetable(iprof, :) = lshape(1)
      ENDIF

      ! layer 1
      IF(nshapelayers .gt. 1) THEN
        lshapetable(iprof, 1:lshape_bot(1)) = lshape(1)
      ENDIF

      !WRITE(*, *) "[Main]: vertinho 2~n-1"
      
      ! layer 2~n-1
      IF(nshapelayers .gt. 2) THEN
        DO ishapelayer = 2, nshapelayers-1
            lshapetable(iprof, lshape_bot(ishapelayer-1)+1:lshape_bot(ishapelayer)) &
            & = lshape(ishapelayer)
        ENDDO
      ENDIF

      !WRITE(*, *) "[Main]: vertinho n"

      ! layer n
      IF(nshapelayers .gt. 1) THEN
        lshapetable(iprof, lshape_bot(nshapelayers-1)+1:nlevels) = lshape(nshapelayers)
      ENDIF

    ENDDO

    !WRITE(*,*) "profile1:"
    !WRITE(*,*) npadlevels(1)
    !WRITE(*,*) lshapetable(1, :)
    !WRITE(*,*) "profile5100:"
    !WRITE(*,*) npadlevels(7673)
    !WRITE(*,*) lshapetable(7673, :)
    !WRITE(*,*) "profile10201:"
    !WRITE(*,*) npadlevels(10201)
    !WRITE(*,*) lshapetable(10201, :)


  ENDIF

  !---------------[B] surface field--------------------

  b_rects = 859 !t2m
  varname = "t2m" 
  call read_surface_field(iomd,model_filename, varname, b_rects, level_raw)
  call select_unstructured_values(level_raw, level_selected_raw, select_table, nprof)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % t  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  b_rects = 858 !q2m
  varname = "q2m"
  call read_surface_field(iomd, model_filename, varname, b_rects, level_raw)
  call select_unstructured_values(level_raw, level_selected_raw, select_table, nprof)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % q  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)
 
  b_rects = 860 !u10m
  varname = "u10m"
  call read_surface_field(iomd, model_filename, varname, b_rects, level_raw)
  call select_unstructured_values(level_raw, level_selected_raw, select_table, nprof)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % u  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  b_rects = 861 !v10m
  varname = "v10m"
  call read_surface_field(iomd, model_filename, varname, b_rects, level_raw)
  call select_unstructured_values(level_raw, level_selected_raw, select_table, nprof)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % v  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

  b_rects = 853 !ts
  varname = "ts"
  call read_surface_field(iomd, model_filename, varname, b_rects, level_raw)
  call select_unstructured_values(level_raw, level_selected_raw, select_table, nprof)
  DO iprof = 1, nprof
        profiles(iprof) % skin % t  = level_selected_raw(iprof)
  ENDDO
  deallocate(level_selected_raw)
  deallocate(level_raw)

!---------------[C] original full level -----------------------

  ! on original full level

  b_rects = 699 ! cldfra
  varname = "cldfra"
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)
  DO iprof=1, nprof
    cld_profiles(iprof)%cc(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  b_rects = 349 ! Qc
  varname = "Qc"
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)
  DO iprof=1, nprof
    cld_profiles(iprof)%clw(:) = levels_selected_raw(:, iprof) 
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  b_rects = 449 ! Qi
  varname = "Qi"
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)
  DO iprof=1, nprof
    cld_profiles(iprof)%ciw(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  b_rects = 399 ! Qr
  varname = "Qr"
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)
  DO iprof=1, nprof
    cld_profiles(iprof)%rain(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  b_rects = 499 ! Qs
  varname = "Qs"
  call read_levels_field(iomd, model_filename, varname, b_rects, levels_raw, nlevels)
  call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nprof)
  DO iprof=1, nprof
    cld_profiles(iprof)%sp(:) = levels_selected_raw(:, iprof)
  ENDDO
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  !---------------[D]  internally generated--------------

  DO iprof = 1, nprof
        profiles(iprof) % skin % salinity  = 35.0
        profiles(iprof) % skin % fastem(:) = (/3.0, 5.0, 15.0, 0.1, 0.3/)
        profiles(iprof) % skin % surftype  = 1
        profiles(iprof) % skin % watertype = 1
        profiles(iprof) % elevation = 0
  ENDDO

  deallocate(select_table)


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


  !---------------------------------------------------------------------------
  !  *. output the profile for a test
  !---------------------------------------------------------------------------

  ! deleted block for profile output_test

  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV-SCATT forward model
  ! --------------------------------------------------------------------------
  IF (nthreads <= 1) THEN
    WRITE(*,*) "before enter rttov_scatt_zvertinho"
    CALL rttov_scatt_zvertinho ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for Mietable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coefs,               &! in    coefficients structure
        nmietables,          &! in    number of mietables
        coefs_scatt_sets,    &! in    Mietable structure(array:(nmietables))
        lshapetable,         &! in    shape index(array:(nprof, nlevels))
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        radiance,            &! inout computed radiances
        gh = gh)              ! in    geopotential height [km]
    WRITE(*,*) "after enter rttov_scatt_zvertinho"
  ELSE
    CALL rttov_parallel_scatt ( &
        errorstatus,         &! out   error flag
        opts_scatt,          &! in    RTTOV-SCATT options structure
        nlevels,             &! in    number of profile levels
        chanprof,            &! in    channel and profile index structure
        frequencies,         &! in    channel indexes for Mietable lookup
        profiles,            &! in    profile array
        cld_profiles,        &! in    cloud/hydrometeor profile array
        coefs,               &! in    coefficients structure
        coefs_scatt_sets(1), &! in    Mietable structure(array:(nmietables))
        calcemis,            &! in    flag for internal emissivity calcs
        emissivity,          &! inout input/output emissivities per channel
        radiance,            &! inout computed radiances
        nthreads = nthreads)  ! in    number of threads to use
  ENDIF
  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_scatt error'
    CALL rttov_exit(errorstatus)
  ENDIF

  !---------------------------------------------------------------------------
  !  *. output the profile for a test
  !---------------------------------------------------------------------------

  ! deleted block for test_output

  !=====================================================
  !============== Output results == start ==============

  ! Open output file where results are written
  OPEN(ioout, file=TRIM(output_filename), status='unknown', form='formatted', iostat=ios)
  IF (ios /= 0) THEN
    WRITE(*,*) 'error opening the output file ios= ', ios
    CALL rttov_exit(errorstatus_fatal)
  ENDIF

  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    WRITE(ioout,333) profiles(iprof)%longitude, profiles(iprof)%latitude
    WRITE(ioout,444) (radiance % bt(j), j = 1+joff, nchannels+joff)

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

  DEALLOCATE(cld_profiles, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'dellocation error for cld_profiles array'
  ENDIF

  ! Deallocate structures for RTTOV direct model
  CALL rttov_alloc_direct( &
        errorstatus,           &
        0_jpim,                &  ! 0 => deallocate
        nprof,                 &
        nchanprof,             &
        nlevels,               &
        chanprof,              &
        opts,                  &
        profiles,              &
        coefs,                 &
        radiance = radiance,   &
        calcemis = calcemis,   &
        emissivity = emissivity)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_direct structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  DO imietable = 1, nmietables
    CALL rttov_dealloc_scattcoeffs(coefs_scatt_sets(imietable))
  ENDDO
  DEALLOCATE(coefs_scatt_sets)

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
222  FORMAT(1X,10F10.5)
333  FORMAT(1X,2F8.3)
444  FORMAT(1X,10F8.3)
777  FORMAT(/,A,A9,I3)

END PROGRAM rttovscatt_fwd_zvertinho_interp_main
