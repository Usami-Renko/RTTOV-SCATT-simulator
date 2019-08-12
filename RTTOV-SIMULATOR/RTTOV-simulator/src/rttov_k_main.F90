PROGRAM rttov_k_main

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
         rttov_coefs,         &
         rttov_profile,       &
         rttov_transmission,  &
         rttov_radiance,      &
         rttov_chanprof,      &
         rttov_emissivity,    &
         rttov_reflectance

  ! jpim, jprb and jplm are the RTTOV integer, real and logical KINDs
  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_unix_env, ONLY : rttov_exit

  ! my mod
  USE rttovscatt_mod, ONLY : plevels, delta

  IMPLICIT NONE

#include "rttov_k.interface"
#include "rttov_parallel_k.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"
#include "rttov_alloc_k.interface"
#include "rttov_init_prof.interface"
#include "rttov_init_rad.interface"
#include "rttov_init_transmission.interface"
#include "rttov_user_options_checkinput.interface"
#include "rttov_print_opts.interface"
#include "rttov_print_profile.interface"
#include "rttov_skipcommentline.interface"

! my include
#include "rttovscatt_sva.interface"
#include "rttovscatt_select_level.interface"
#include "rttovscatt_select_levels.interface"
#include "rttovscatt_gribapi.interface"
#include "rttovscatt_comput.interface"

  !--------------------------
  !
  INTEGER(KIND=jpim), PARAMETER :: ioin   = 20        ! unit for input sva file
  INTEGER(KIND=jpim), PARAMETER :: ioout  = 21        ! unit for output
  INTEGER(KIND=jpim), PARAMETER :: ioouttest  = 22    ! unit for output test profile

  ! RTTOV variables/structures
  !====================
  TYPE(rttov_options)              :: opts                       ! Options structure
  TYPE(rttov_coefs)                :: coefs                      ! Coefficients structure
  TYPE(rttov_chanprof),    POINTER :: chanprof(:)      => NULL() ! Input channel/profile list
  LOGICAL(KIND=jplm),      POINTER :: calcemis(:)      => NULL() ! Flag to indicate calculation of emissivity within RTTOV
  TYPE(rttov_emissivity),  POINTER :: emissivity(:)    => NULL() ! Input/output surface emissivity
  TYPE(rttov_emissivity),  POINTER :: emissivity_k(:)  => NULL() ! Emissivity Jacobians
  LOGICAL(KIND=jplm),      POINTER :: calcrefl(:)      => NULL() ! Flag to indicate calculation of BRDF within RTTOV
  TYPE(rttov_reflectance), POINTER :: reflectance(:)   => NULL() ! Input/output surface BRDF
  TYPE(rttov_reflectance), POINTER :: reflectance_k(:) => NULL() ! Reflectance Jacobians
  TYPE(rttov_profile),     POINTER :: profiles(:)      => NULL() ! Input profiles
  TYPE(rttov_profile),     POINTER :: profiles_k(:)    => NULL() ! Output Jacobians
  TYPE(rttov_transmission)         :: transmission               ! Output transmittances
  TYPE(rttov_transmission)         :: transmission_k             ! Transmittance Jacobians
  TYPE(rttov_radiance)             :: radiance                   ! Output radiances
  TYPE(rttov_radiance)             :: radiance_k                 ! Radiance Jacobians

  INTEGER(KIND=jpim)               :: errorstatus                ! Return error status of RTTOV subroutine calls

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

  CHARACTER(LEN=32)                             :: shortName

  ! variables for input
  !====================
  CHARACTER(LEN=256) :: coef_filename
  CHARACTER(LEN=256) :: sva_filename
  CHARACTER(LEN=256) :: model_filename
  CHARACTER(LEN=256) :: output_filename
  INTEGER(KIND=jpim) :: nthreads
  INTEGER(KIND=jpim) :: dosolar
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: nprof
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nchanprof
  INTEGER(KIND=jpim), ALLOCATABLE :: channel_list(:)
  REAL(KIND=jprb)    :: trans_out(10)
  CHARACTER(LEN=11)  :: gas_unit
  ! loop variables
  INTEGER(KIND=jpim) :: j, jch, l
  INTEGER(KIND=jpim) :: np, nch
  INTEGER(KIND=jpim) :: ilev, nprint
  INTEGER(KIND=jpim) :: iprof, joff
  INTEGER            :: ios

  !- End of header --------------------------------------------------------

  ! The usual steps to take when running RTTOV are as follows:
  !   1. Specify required RTTOV options
  !   2. Read coefficients
  !   3. Allocate RTTOV input and output structures
  !   4. Set up the chanprof array with the channels/profiles to simulate
  !   5. Read input profile(s)
  !   6. Set up surface emissivity and/or reflectance
  !   7. Call rttov_k and store results
  !   8. Deallocate all structures and arrays

  ! If nthreads is greater than 1 the parallel RTTOV interface is used.
  ! To take advantage of multi-threaded execution you must have compiled
  ! RTTOV with openmp enabled. See the user guide and the compiler flags.

  errorstatus = 0_jpim

  !=====================================================
  !========== Interactive inputs == start ==============

  WRITE(0,*) 'enter path of coefficient file'
  READ(*,*) coef_filename
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
  WRITE(0,*) 'turn on solar simulations? (0=no, 1=yes)'
  READ(*,*) dosolar
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


  ! --------------------------------------------------------------------------
  ! 1. Initialise RTTOV options structure
  ! --------------------------------------------------------------------------

  IF (dosolar == 1) THEN
    opts % rt_ir % addsolar = .TRUE.           ! Include solar radiation
  ELSE
    opts % rt_ir % addsolar = .FALSE.          ! Do not include solar radiation
  ENDIF
  opts % interpolation % addinterp   = .TRUE.  ! Allow interpolation of input profile
  opts % interpolation % interp_mode = 1       ! Set interpolation method
  opts % rt_all % addrefrac          = .TRUE.  ! Include refraction in path calc
  opts % rt_ir % addclouds           = .FALSE. ! Don't include cloud effects
  opts % rt_ir % addaerosl           = .FALSE. ! Don't include aerosol effects

  opts % rt_ir % ozone_data          = .FALSE. ! Set the relevant flag to .TRUE.
  opts % rt_ir % co2_data            = .FALSE. !   when supplying a profile of the
  opts % rt_ir % n2o_data            = .FALSE. !   given trace gas (ensure the
  opts % rt_ir % ch4_data            = .FALSE. !   coef file supports the gas)
  opts % rt_ir % co_data             = .FALSE. !
  opts % rt_ir % so2_data            = .FALSE. !
  opts % rt_mw % clw_data            = .FALSE. !

!   opts % rt_all % switchrad          = .FALSE. ! Input K perturbation in radiance
  opts % rt_all % switchrad          = .TRUE.  ! Input K perturbation in BT

  opts % config % verbose            = .TRUE.  ! Enable printing of warnings

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

  ! Ensure the options and coefficients are consistent
  CALL rttov_user_options_checkinput(errorstatus, opts, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'error in rttov options'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 3. Allocate RTTOV input and output structures
  ! --------------------------------------------------------------------------

  ! Determine the total number of radiances to simulate (nchanprof).
  ! In this example we simulate all specified channels for each profile, but
  ! in general one can simulate a different number of channels for each profile.

  nchanprof = nchannels * nprof

  ! Allocate structures for rttov_k
  CALL rttov_alloc_k( &
        errorstatus,                 &
        1_jpim,                      &  ! 1 => allocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        profiles_k,                  &
        coefs,                       &
        transmission,                &
        transmission_k,              &
        radiance,                    &
        radiance_k,                  &
        calcemis=calcemis,           &
        emissivity=emissivity,       &
        emissivity_k=emissivity_k,   &
        calcrefl=calcrefl,           &
        reflectance=reflectance,     &
        reflectance_k=reflectance_k, &
        init=.TRUE._jplm)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'allocation error for rttov_k structures'
    CALL rttov_exit(errorstatus)
  ENDIF


  ! --------------------------------------------------------------------------
  ! 4. Build the list of profile/channel indices in chanprof
  ! --------------------------------------------------------------------------

  nch = 0_jpim
  DO j = 1, nprof
    DO jch = 1, nchannels
      nch = nch + 1_jpim
      chanprof(nch)%prof = j
      chanprof(nch)%chan = channel_list(jch)
    ENDDO
  ENDDO


  ! --------------------------------------------------------------------------
  ! 5. Read profile data
  ! --------------------------------------------------------------------------

  !===============================================
  !========== Read profiles == start =============
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


  !---------------[B] surface field--------------------

  shortName = "sp"
  call read_surface_field_latlon(model_filename, shortName, 0, level_raw, lats_raw, lons_raw)
  level_raw = level_raw / 100.0
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  call select_values(lats_raw, lats_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  call select_values(lons_raw, lons_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  DO iprof = 1, nprof
        profiles(iprof) % s2m % p  = level_selected_raw(iprof)
        profiles(iprof) % latitude  =  lats_selected_raw(iprof) 
        profiles(iprof) % longitude =  lons_selected_raw(iprof)
  ENDDO
  deallocate(level_raw)
  deallocate(level_selected_raw)
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

  DO iprof = 1, nprof
      profiles(iprof) % p(:) = plevels
  ENDDO

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


  !========== Read profiles == end =============
  !=============================================


  ! --------------------------------------------------------------------------
  ! 6. Specify surface emissivity and reflectance
  ! --------------------------------------------------------------------------

  ! In this example we have no values for input emissivities
  emissivity(:) % emis_in = 0._jprb

  ! Calculate emissivity within RTTOV where the input emissivity value is
  ! zero or less (all channels in this case)
  calcemis(:) = (emissivity(:) % emis_in <= 0._jprb)

  ! In this example we have no values for input reflectances
  reflectance(:) % refl_in = 0._jprb

  ! Calculate BRDF within RTTOV where the input BRDF value is zero or less
  ! (all channels in this case)
  calcrefl(:) = (reflectance(:) % refl_in <= 0._jprb)

  ! Use default cloud top BRDF for simple cloud in VIS/NIR channels
  reflectance(:) % refl_cloud_top = 0._jprb


  ! --------------------------------------------------------------------------
  ! 7. Call RTTOV K model
  ! --------------------------------------------------------------------------

  ! The input/output K variables must be initialised to zero before every call to rttov_k:

  ! Initialise profiles_k, radiance_k and transmission_k structures to zero
  CALL rttov_init_prof(profiles_k(:))
  CALL rttov_init_rad(radiance_k)
  CALL rttov_init_transmission(transmission_k)

  ! Initialise emissivity_k to zero
  emissivity_k(:) % emis_in = 0._jprb         ! Output emissivity Jacobian is in emis_in
  emissivity_k(:) % emis_out = 0._jprb

  ! Initialise reflectance_k to zero
  reflectance_k(:) % refl_in = 0._jprb        ! Output reflectance Jacobian is in refl_in
  reflectance_k(:) % refl_out = 0._jprb
  reflectance_k(:) % refl_cloud_top = 0._jprb ! This is not an active variable in the TL/AD/K

  ! Set input perturbation in radiance_k:
  !   If switchrad is TRUE the perturbation in bt(:) is used for "thermal" channels
  !   If switchrad is FALSE the perturbation in total(:) is used
  !   For solar-only channels the perturbation in total(:) is *always* used
  !   It is harmless to specify inputs in both bt/total in any case: RTTOV will use the
  !     appropriate input for each channel
  radiance_k % total(:) = 1._jprb
  radiance_k % bt(:) = 1._jprb

  IF (nthreads <= 1) THEN
    CALL rttov_k(                          &
            errorstatus,                   &! out   error flag
            chanprof,                      &! in    channel and profile index structure
            opts,                          &! in    options structure
            profiles,                      &! in    profile array
            profiles_k,                    &! inout Jacobian array
            coefs,                         &! in    coefficients structure
            transmission,                  &! inout computed transmittances
            transmission_k,                &! inout transmittance Jacobians
            radiance,                      &! inout computed radiances
            radiance_k,                    &! inout input radiance/BT perturbation
            calcemis      = calcemis,      &! in    flag for internal emissivity calcs
            emissivity    = emissivity,    &! inout input/output emissivities per channel
            emissivity_k  = emissivity_k,  &! inout emissivity Jacobians
            calcrefl      = calcrefl,      &! in    flag for internal BRDF calcs
            reflectance   = reflectance,   &! inout input/output BRDFs per channel
            reflectance_k = reflectance_k)  ! inout BRDF Jacobians
  ELSE
    CALL rttov_parallel_k(                 &
            errorstatus,                   &! out   error flag
            chanprof,                      &! in    channel and profile index structure
            opts,                          &! in    options structure
            profiles,                      &! in    profile array
            profiles_k,                    &! inout Jacobian array
            coefs,                         &! in    coefficients structure
            transmission,                  &! inout computed transmittances
            transmission_k,                &! inout transmittance Jacobians
            radiance,                      &! inout computed radiances
            radiance_k,                    &! inout input radiance/BT perturbation
            calcemis      = calcemis,      &! in    flag for internal emissivity calcs
            emissivity    = emissivity,    &! inout input/output emissivities per channel
            emissivity_k  = emissivity_k,  &! inout emissivity Jacobians
            calcrefl      = calcrefl,      &! in    flag for internal BRDF calcs
            reflectance   = reflectance,   &! inout input/output BRDFs per channel
            reflectance_k = reflectance_k, &! inout BRDF Jacobians
            nthreads      = nthreads)       ! in    number of threads to use
  ENDIF

  IF (errorstatus /= errorstatus_success) THEN
    WRITE (*,*) 'rttov_k error'
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

  CALL rttov_print_profile(profiles(1), lu=ioouttest)
  CALL rttov_print_profile(profiles(121), lu=ioouttest)

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
    gas_unit = '(K/(kg/kg))'
  ELSE
    gas_unit = '(K/ppmv)   '
  ENDIF


  DO iprof = 1, nprof

    joff = (iprof-1_jpim) * nchannels

    ! WRITE(ioout,*)' Profile ', iprof

    ! CALL rttov_print_profile(profiles(iprof), lu=ioout)
    ! CALL rttov_print_cld_profile(cld_profiles(iprof), lu=ioout)

    WRITE(ioout,444) (radiance % bt(j), j = 1+joff, nchannels+joff)


    DO j = 1, nchannels

      WRITE(ioout,*)'Channel ', chanprof(j+joff) % chan
      WRITE(ioout,'(/a5,a9,a18,1x,a19)',advance='no') "Level", "Pressure", "JAC(Temp) (K/K)", "JAC(WV) "//gas_unit
      IF (opts % rt_ir % ozone_data) THEN
        WRITE(ioout,'(1x,a19)',advance='no') "JAC(O3) "//gas_unit
      ENDIF
      IF (opts % rt_ir % co2_data) THEN
        WRITE(ioout,'(1x,a20)',advance='no') "JAC(CO2) "//gas_unit
      ENDIF
      IF (opts % rt_ir % n2o_data) THEN
        WRITE(ioout,'(1x,a20)',advance='no') "JAC(N2O) "//gas_unit
      ENDIF
      IF (opts % rt_ir % co_data) THEN
        WRITE(ioout,'(1x,a19)',advance='no') "JAC(CO) "//gas_unit
      ENDIF
      IF (opts % rt_ir % ch4_data) THEN
        WRITE(ioout,'(1x,a20)',advance='no') "JAC(CH4) "//gas_unit
      ENDIF
      WRITE(ioout,'(a)',advance='yes')

      DO l = 1, profiles_k(j+joff) % nlevels
        WRITE(ioout,'(i4,1x,f9.4,1x,e16.6,1x,e17.6)',advance='no') &
          l, profiles(iprof)%p(l), profiles_k(j+joff)%t(l), profiles_k(j+joff)%q(l)
        IF (opts % rt_ir % ozone_data) THEN
          WRITE(ioout,'(1x,e18.6)',advance='no') profiles_k(j+joff)%o3(l)
        ENDIF
        IF (opts % rt_ir % co2_data) THEN
          WRITE(ioout,'(1x,e19.6)',advance='no') profiles_k(j+joff)%co2(l)
        ENDIF
        IF (opts % rt_ir % n2o_data) THEN
          WRITE(ioout,'(1x,e19.6)',advance='no') profiles_k(j+joff)%n2o(l)
        ENDIF
        IF (opts % rt_ir % co_data) THEN
          WRITE(ioout,'(1x,e18.6)',advance='no') profiles_k(j+joff)%co(l)
        ENDIF
        IF (opts % rt_ir % ch4_data) THEN
          WRITE(ioout,'(1x,e19.6)',advance='no') profiles_k(j+joff)%ch4(l)
        ENDIF
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
  DEALLOCATE (channel_list, stat=alloc_status)
  IF (alloc_status /= 0) THEN
    WRITE(*,*) 'mem dellocation error'
  ENDIF

  ! Deallocate structures for rttov_k
  CALL rttov_alloc_k( &
        errorstatus,                 &
        0_jpim,                      &  ! 0 => deallocate
        nprof,                       &
        nchanprof,                   &
        nlevels,                     &
        chanprof,                    &
        opts,                        &
        profiles,                    &
        profiles_k,                  &
        coefs,                       &
        transmission,                &
        transmission_k,              &
        radiance,                    &
        radiance_k,                  &
        calcemis=calcemis,           &
        emissivity=emissivity,       &
        emissivity_k=emissivity_k,   &
        calcrefl=calcrefl,           &
        reflectance=reflectance,     &
        reflectance_k=reflectance_k)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'deallocation error for rttov_k structures'
    CALL rttov_exit(errorstatus)
  ENDIF

  CALL rttov_dealloc_coefs(errorstatus, coefs)
  IF (errorstatus /= errorstatus_success) THEN
    WRITE(*,*) 'coefs deallocation error'
  ENDIF


! Format definitions for output
111  FORMAT(1X,10I8)
1115 FORMAT(3X,10I8)
222  FORMAT(1X,10F8.2)
444  FORMAT(1X,10F8.3)
4444 FORMAT(1X,10F8.4)
4445 FORMAT(1X,I2,10F8.4)
777  FORMAT(/,A,A9,I3)

END PROGRAM rttov_k_main
