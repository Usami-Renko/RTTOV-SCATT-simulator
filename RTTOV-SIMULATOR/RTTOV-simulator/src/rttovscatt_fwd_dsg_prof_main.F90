PROGRAM rttovscatt_fwd_dsg_prof_main

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
  
    ! my mod
    USE rttovscatt_mod, ONLY : plevels
  
    IMPLICIT NONE
  
#include "rttov_scatt.interface"
#include "rttov_scatt_vertinho_out.interface"
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
#include "rttovscatt_gribapi.interface"
#include "rttovscatt_comput.interface"
#include "rttov_print_opts.interface"
  
  
    !--------------------------
    !
    INTEGER(KIND=jpim), PARAMETER :: ioin   = 20    ! unit for input sva file
    INTEGER(KIND=jpim), PARAMETER :: ioout  = 21    ! unit for output
  
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
    real(KIND=jprb), allocatable    :: levels_selected_raw(:, :)  !(nlevels, nprof)
    real(KIND=jprb), allocatable    :: level_raw(:)                 !(npoints)
    real(KIND=jprb), allocatable    :: level_selected_raw(:)        !(nprof)
  
    ! interp module modified
    real(KIND=jprb), allocatable    :: sva_values(:, :)
    integer(KIND=jpim), allocatable :: select_table(:, :)
    
    real(KIND=jprb), dimension(:), allocatable    :: lats_raw ! (npoints)
    real(KIND=jprb), dimension(:), allocatable    :: lons_raw ! (npoints)
    real(KIND=jprb), dimension(:), allocatable    :: lats_selected_raw ! (npoints)
    real(KIND=jprb), dimension(:), allocatable    :: lons_selected_raw ! (npoints)
  
    real(KIND=jprb), dimension(:), allocatable    :: ph
    CHARACTER(LEN=32)                             :: shortName
  
    ! padding variables
    integer(KIND=jpim)                            :: ipadlevels
    integer(KIND=jpim)                            :: npad
    real(KIND=jprb), dimension(:), allocatable    :: single_ph
    real(KIND=jprb)                               :: pad_increment = 0.01
    real(KIND=jprb)                               :: topT
    real(KIND=jprb)                               :: topq
  
  
    ! variables for input
    !====================
    CHARACTER(LEN=256) :: coef_filename
    CHARACTER(LEN=256) :: mietable_dir
    CHARACTER(LEN=256) :: mietable_filename
    CHARACTER(LEN=256) :: sva_filename
    CHARACTER(LEN=256) :: model_filename
    CHARACTER(LEN=256) :: output_dir
    CHARACTER(LEN=256) :: avgprof_filename
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
  
    ! variables for input designed profile
    INTEGER(KIND=jpim) :: H_ngrid, L_ngrid, H_igrid, L_igrid, ilevel 
    REAL(KIND=jprb), dimension(:), allocatable  :: H_grid, L_grid
    REAL(KIND=jprb), dimension(:), allocatable  :: avgprof, tmpavgprof
    REAL(KIND=jprb), dimension(:,:,:), allocatable  :: packed_out
    REAL(KIND=jprb)    :: mycfrac 

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
    ! WRITE(0,*) 'enter path of OUTPUT directory'
    READ(*,*) output_dir
    ! WRITE(0,*) 'enter path of AvgProf file'
    READ(*,*) avgprof_filename
    ! WRITE(0,*) 'enter number of profile levels'
    READ(*,*) nlevels
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

    ! Input variables of designed profile
    ! WRITE(0,*) 'enter designed High ICE cloud ngrid'
    READ(*,*) H_ngrid
    ALLOCATE(H_grid(H_ngrid))
    ! WRITE(0,*) 'enter designed High ICE cloud grid'
    READ(*,*,iostat=ios) H_grid
    ! WRITE(0,*) 'enter designed Low  ICE cloud ngrid'
    READ(*,*) L_ngrid
    ALLOCATE(L_grid(L_ngrid))
    ! WRITE(0,*) 'enter designed Low ICE cloud grid'
    READ(*,*,iostat=ios) L_grid
    nprof = H_ngrid * L_ngrid
  
    ! input variables of vertical inhomogeneity
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
    opts_scatt % lusercfrac = .True.
  
  
  
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
          nlevels,                 &
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
  
    !--------------[A0] SVA---------------------------------
  
    call read_unstructured_sva_fields(ioin, sva_filename, sva_values, select_table, 1)
  
    DO iprof = 1, nprof
        profiles(iprof) % zenangle = sva_values(1, 1)
        profiles(iprof) % azangle = sva_values(2, 1)
    ENDDO
  
    deallocate(sva_values)
  
    !---------------[A1] p & ph & populate npadlevels & sp--------------------
  
    allocate(ph(nlevels))         ! ph(1:nlevels)
    allocate(single_ph(nlevels))  ! ph(1:nlevels)
  
    call generate_std_levels(plevels, ph)
  
    shortName = "sp"
    call read_surface_field_latlon(model_filename, shortName, 0, level_raw, lats_raw, lons_raw)
    level_raw = level_raw / 100.0
    call select_unstructured_values(level_raw, level_selected_raw, select_table, 1)
    call select_unstructured_values(lats_raw, lats_selected_raw, select_table, 1)
    call select_unstructured_values(lons_raw, lons_selected_raw, select_table, 1)
  
    !* get npad for the single profile
    DO ipadlevels = 0, nlevels-1
      IF (ph(nlevels-ipadlevels) .lt. level_selected_raw(1)) THEN
        npad = ipadlevels
        EXIT
      ENDIF
    ENDDO

    DO iprof = 1, nprof
        !---------- populate profiles(:)%p(npadlevels+1:nlevels)
        profiles(iprof)%p(npad+1:nlevels) = plevels(1:nlevels-npad)
        !---------- populate padding levels : profiles(:)%p(1:npadlevels) 
        DO ipadlevels = 1, npad
          profiles(iprof)%p(ipadlevels) = 0.005 + ipadlevels*pad_increment
        ENDDO
        !---------- populate profiles(:)%ph(1:nlevels)
        IF (iprof .eq. 1) THEN 
          CALL generate_std_levels(profiles(iprof)%p, single_ph)
        ENDIF
        cld_profiles(iprof)%ph(1:nlevels) = single_ph(:) 
        !-----------other variables
        profiles(iprof) % s2m % p  = level_selected_raw(1)
        cld_profiles(iprof) % ph(nlevels+1) = level_selected_raw(1)
        profiles(iprof) % latitude  =  lats_selected_raw(1) 
        profiles(iprof) % longitude =  lons_selected_raw(1)
    ENDDO
  
    deallocate(ph)
    deallocate(level_selected_raw)
    deallocate(level_raw)
    deallocate(lats_raw)
    deallocate(lons_raw)
    deallocate(lats_selected_raw)
    deallocate(lons_selected_raw)
  
  
    !---------------[A2] populate the lshapetable---------
    ! vertical inhomogeniety mode 1 : shapelayers
  
    IF(vertinho_mode .eq. 1) THEN
  
      allocate(lshapetable(nprof, nlevels))
  
      DO iprof = 1, nprof
  
        ! IF (npad .gt. 1) THEN
        !   WRITE(*,*) npad, iprof
        ! ENDIF
  
        !WRITE(*, *) "[Main]: vertinho 1"
  
        IF(nshapelayers .eq. 1) THEN
          lshapetable(iprof, :) = lshape(1)
        ENDIF
  
        ! layer 1
        IF(nshapelayers .gt. 1) THEN
          lshapetable(iprof, 1:lshape_bot(1)+npad) = lshape(1)
        ENDIF
  
        !WRITE(*, *) "[Main]: vertinho 2~n-1"
        
        ! layer 2~n-1
        IF(nshapelayers .gt. 2) THEN
          DO ishapelayer = 2, nshapelayers-1
              lshapetable(iprof, lshape_bot(ishapelayer-1)+npad+1:lshape_bot(ishapelayer)+npad) &
              & = lshape(ishapelayer)
          ENDDO
        ENDIF
  
        !WRITE(*, *) "[Main]: vertinho n"
  
        ! layer n
        IF(nshapelayers .gt. 1) THEN
          lshapetable(iprof, lshape_bot(nshapelayers-1)+npad+1:nlevels) = lshape(nshapelayers)
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
  
    shortName = "2t"
    call read_surface_field(model_filename, shortName, 2, level_raw)
    call select_unstructured_values(level_raw, level_selected_raw, select_table, 1)
    DO iprof = 1, nprof
          profiles(iprof) % s2m % t  = level_selected_raw(1)
    ENDDO
    deallocate(level_selected_raw)
    deallocate(level_raw)
  
    shortName = "q"
    call read_surface_field(model_filename, shortName, 2, level_raw)
    call select_unstructured_values(level_raw, level_selected_raw, select_table, 1)
    DO iprof = 1, nprof
          profiles(iprof) % s2m % q  = level_selected_raw(1)
    ENDDO
    deallocate(level_selected_raw)
    deallocate(level_raw)
  
    shortName = "10u"
    call read_surface_field(model_filename, shortName, 10, level_raw)
    call select_unstructured_values(level_raw, level_selected_raw, select_table, 1)
    DO iprof = 1, nprof
          profiles(iprof) % s2m % u  = level_selected_raw(1)
    ENDDO
    deallocate(level_selected_raw)
    deallocate(level_raw)
  
    shortName = "10v"
    call read_surface_field(model_filename, shortName, 10, level_raw)
    call select_unstructured_values(level_raw, level_selected_raw, select_table, 1)
    DO iprof = 1, nprof
          profiles(iprof) % s2m % v  = level_selected_raw(1)
    ENDDO
    deallocate(level_selected_raw)
    deallocate(level_raw)
  
    shortName = "t"
    call read_surface_field(model_filename, shortName, 0, level_raw)
    call select_unstructured_values(level_raw, level_selected_raw, select_table, 1)
    DO iprof = 1, nprof
          profiles(iprof) % skin % t  = level_selected_raw(1)
    ENDDO
    deallocate(level_selected_raw)
    deallocate(level_raw)
  
  !---------------[C]  levels field-----------------------
  
    WRITE(*,*) 'Before levels field'
  
    shortName = "t"
    call read_levels_field(model_filename, shortName, levels_raw, nlevels)
    call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, 1)
    DO iprof=1, nprof
      topT = levels_selected_raw(1, 1)
      profiles(iprof)%t(npad+1:nlevels) = levels_selected_raw(1:nlevels-npad, 1)
      profiles(iprof)%t(1:npad) = topT
    ENDDO
    deallocate(levels_selected_raw)
    deallocate(levels_raw)
  
    shortName = "q"
    call read_levels_field(model_filename, shortName, levels_raw, nlevels)
    call select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, 1)
    DO iprof=1, nprof
      topq = levels_selected_raw(1, 1)
      profiles(iprof)%q(npad+1:nlevels) = levels_selected_raw(1:nlevels-npad, 1)
      profiles(iprof)%q(1:npad) = topq
    ENDDO
    deallocate(levels_selected_raw)
    deallocate(levels_raw)

    ! write(*,*) profiles(1)%q
    ! write(*,*) profiles(1)%t
    ! write(*,*) profiles(1)%s2m%u
    ! write(*,*) profiles(1)%s2m%v
    ! write(*,*) profiles(1)%skin%t
    ! write(*,*) profiles(1)%s2m%t
    ! write(*,*) profiles(1)%s2m%q

    !---------------[D]  read from avgprof.dat -----------------------
    allocate(avgprof(nlevels))
    allocate(tmpavgprof(nlevels))
    open(ioin, file=trim(avgprof_filename), status='old', iostat=ios)

    if (ios /= 0) then
    	WRITE(*,*) 'error opening avgprof file ios= ', ios
    endif
    
    !* cc
    read(ioin, *) avgprof
    DO H_igrid = 1, H_ngrid
      DO L_igrid = 1, L_ngrid
        iprof = (H_igrid - 1) * L_ngrid + L_igrid
        cld_profiles(iprof)%cc(npad+1:nlevels) = avgprof(1:nlevels-npad) / 100.
        cld_profiles(iprof)%cc(1:npad) = 0.0
      ENDDO
    ENDDO

    !* ciw
    read(ioin, *) avgprof
    DO H_igrid = 1, H_ngrid
      tmpavgprof(1:14) = avgprof(1:14) * H_grid(H_igrid) 
      DO L_igrid = 1, L_ngrid
        tmpavgprof(15:30) = avgprof(15:30) * L_grid(L_igrid)

        iprof = (H_igrid - 1) * L_ngrid + L_igrid
        cld_profiles(iprof)%ciw(npad+1:nlevels) = tmpavgprof(1:nlevels-npad)
        cld_profiles(iprof)%ciw(1:npad) = 0.0
      ENDDO
    ENDDO

    !* clw
    read(ioin, *) avgprof
    DO H_igrid = 1, H_ngrid
      DO L_igrid = 1, L_ngrid  
        iprof = (H_igrid - 1) * L_ngrid + L_igrid
        cld_profiles(iprof)%clw(npad+1:nlevels) = avgprof(1:nlevels-npad)
        cld_profiles(iprof)%clw(1:npad) = 0.0
      ENDDO
    ENDDO

    !* rain
    read(ioin, *) avgprof
    DO H_igrid = 1, H_ngrid
      DO L_igrid = 1, L_ngrid  
        iprof = (H_igrid - 1) * L_ngrid + L_igrid
        cld_profiles(iprof)%rain(npad+1:nlevels) = avgprof(1:nlevels-npad)
        cld_profiles(iprof)%rain(1:npad) = 0.0
      ENDDO
    ENDDO

    !* sp
    read(ioin, *) avgprof
    DO H_igrid = 1, H_ngrid
      tmpavgprof(1:14) = avgprof(1:14) * H_grid(H_igrid) 
      DO L_igrid = 1, L_ngrid
        tmpavgprof(15:30) = avgprof(15:30) * L_grid(L_igrid)

        iprof = (H_igrid - 1) * L_ngrid + L_igrid
        cld_profiles(iprof)%sp(npad+1:nlevels) = tmpavgprof(1:nlevels-npad)
        cld_profiles(iprof)%sp(1:npad) = 0.0
      ENDDO
    ENDDO

    read(ioin, *) mycfrac
    DO H_igrid = 1, H_ngrid 
      DO L_igrid = 1, L_ngrid
        iprof = (H_igrid - 1) * L_ngrid + L_igrid
        cld_profiles(iprof)%cfrac = mycfrac / 100.
      ENDDO
    ENDDO

    close(ioin)
    deallocate(avgprof)
  
    !---------------[E]  internally generated--------------
  
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
      WRITE(*,*) "before enter rttov_scatt_vertinho_out"
      CALL rttov_scatt_vertinho_out ( &
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
          packed_out)           ! out   irad_do, irad_up, j_do, j_up, tau (5, nchannels, nlevels+1)    
      WRITE(*,*) "after enter rttov_scatt_vertinho_out"
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
    OPEN(ioout, file=TRIM(output_dir)//"/bt.dat", status='unknown', form='formatted', iostat=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) 'error opening the output file ios= ', ios
      CALL rttov_exit(errorstatus_fatal)
    ENDIF
  
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      WRITE(ioout,444) (radiance % bt(j), j = 1+joff, nchannels+joff)
    ENDDO
  
    ! Close output file
    CLOSE(ioout, iostat=ios)
    IF (ios /= 0) THEN
      WRITE(*,*) 'error closing the output file ios= ', ios
      CALL rttov_exit(errorstatus_fatal)
    ENDIF

    OPEN(ioout, file=TRIM(output_dir)//"/irad_do.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels + 1
        WRITE(ioout,999) (packed_out(1, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/irad_up.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels + 1
        WRITE(ioout,999) (packed_out(2, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/j_do.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels
        WRITE(ioout,999) (packed_out(3, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/j_up.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels
        WRITE(ioout,999) (packed_out(4, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/tau.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels
        WRITE(ioout,999) (packed_out(5, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/ext.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels
        WRITE(ioout,999) (packed_out(6, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/ssa.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels
        WRITE(ioout,999) (packed_out(7, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    OPEN(ioout, file=TRIM(output_dir)//"/asm.dat", status='unknown', form='formatted', iostat=ios)
    DO iprof = 1, nprof
      joff = (iprof-1_jpim) * nchannels
      DO ilevel = 1, nlevels
        WRITE(ioout,999) (packed_out(8, j, ilevel), j = 1+joff, nchannels+joff)
      ENDDO 
    ENDDO
    CLOSE(ioout, iostat=ios)

    deallocate(packed_out)

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
  999  FORMAT(1X,10E15.6)
  
  END PROGRAM rttovscatt_fwd_dsg_prof_main
  