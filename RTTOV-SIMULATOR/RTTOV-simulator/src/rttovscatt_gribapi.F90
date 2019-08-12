! internally generated : FASTEM, Salinity, Surface type, watertype /
!                        Elevation, latitude, longitude

! Varname we use shortName
! we only get "values" in gribapi

!---------------------------------------------------------------

! profile(nprofiles)%[level_var](nlevels)
! profile(nprofiles)%s2m%[s2m_var]
! profile(nprofiles)%skin%[skin_var]

! =====================GRIBAPI_subroutines==============================


! t,      q,      cc,     clw,     ciw,      rain,       sp;
! t,      q,     ccl,   clwmr,    icmr,      rwmr,     snmr;

subroutine read_levels_field(modelfilename, varname, values, nlevels)
  use parkind1, only: jpim, jprb
  use rttovscatt_mod, only: npoints, level_index_start, level_index_stop
  use grib_api

  implicit none

#include "test_MinMax.interface"

  ! dummy arguments
  character(LEN=256), intent(in)                             :: modelfilename
  character(LEN=32), intent(in)                              :: varname
  integer(KIND=jpim), intent(in)                             :: nlevels
  real(KIND=jprb), intent(out), dimension(:, :), allocatable :: values


  ! local variables
  character(LEN=*), parameter        :: open_mode="r"
  integer (KIND=jpim)                :: ilev
  integer (KIND=jpim)                :: iloop
  integer (KIND=jpim)                :: iindex
  integer (KIND=jpim)                :: ifile
  integer (KIND=jpim)                :: istatus
  integer (KIND=jpim), allocatable   :: igrib(:)
  integer (KIND=jpim)                :: level_size
  integer (KIND=jpim), allocatable   :: level_table(:)
  real (KIND=jprb), allocatable      :: lats(:)
  real (KIND=jprb), allocatable      :: lons(:)
  real (KIND=jprb)                   :: mymax_test
  real (KIND=jprb)                   :: mymin_test
  real (KIND=jprb), allocatable      :: values_tmp(:)


  allocate(values(nlevels, npoints))
  allocate(lats(npoints))
  allocate(lons(npoints))
  allocate(values_tmp(npoints))

  call grib_open_file(ifile, modelfilename, open_mode, istatus)

  call grib_index_create(iindex, modelfilename, 'shortName,level')

  call grib_index_get_size(iindex, 'level', level_size)

  allocate(level_table(level_size))
  allocate(igrib(level_size))

  call grib_index_get(iindex, 'level', level_table)

  DO iloop = level_index_start, level_index_stop
    call grib_index_select(iindex, 'shortName', varname)
    call grib_index_select(iindex, 'level', level_table(iloop))
    call grib_new_from_index(iindex, igrib(iloop), istatus)
  ENDDO

  DO iloop = level_index_start, level_index_stop
    ilev = iloop-level_index_start + 1

    ! pass it
    call grib_get_data(igrib(iloop), lats, lons, values_tmp, istatus)
    values(ilev, :) = values_tmp

    ! test
    call mymax(values(ilev,:), npoints, mymax_test)
    call mymin(values(ilev,:), npoints, mymin_test)
    write(*,*) "shortName=", varname, "level=", level_table(iloop), "max=", mymax_test, "min=", mymin_test
  ENDDO


  DO iloop = level_index_start, level_index_stop
    call grib_release(igrib(iloop), istatus)
  ENDDO

  call grib_index_release(iindex, istatus)
  call grib_close_file(ifile, istatus)

  ! deallocate

  deallocate(igrib)
  deallocate(level_table)
  deallocate(lats)
  deallocate(lons)
  deallocate(values_tmp)

end subroutine read_levels_field

!------------------------------------------------------------------------

! we take ps as 2mp
! ph(nlevels+1)(,i.e. 2mp), 2mT, 2mq, 10mU, 10mV, skinT;
! shortName             sp,  2t,   q,  10u,  10v,     t;
! count                393,  12,  11,   13,   14,     5;
! level                  0,   2,   2,   10,   10,     0;
! multi                   ,    ,    ,    *,    *,     *;
! name                   p,   t,    q,   u,    v, skin%t;

subroutine read_surface_field(modelfilename, varname, level, values)
  use parkind1, only : jpim, jprb
  use rttovscatt_mod, only : npoints
  use grib_api

  implicit none

#include "test_MinMax.interface"

  ! dummy arguments
  character(LEN=256), intent(in)                             :: modelfilename
  character(LEN=32), intent(in)                              :: varname
  integer(KIND=jpim), intent(in)                             :: level
  real(KIND=jprb), intent(out), dimension(:), allocatable    :: values ! (npoints)

  ! local variables
  character(LEN=*), parameter        :: open_mode="r"
  integer (KIND=jpim)                :: iindex
  integer (KIND=jpim)                :: ifile
  integer (KIND=jpim)                :: istatus
  integer (KIND=jpim)                :: igrib
  real (KIND=jprb), allocatable      :: lats(:)
  real (KIND=jprb), allocatable      :: lons(:)
  real (KIND=jprb)                   :: mymax_test
  real (KIND=jprb)                   :: mymin_test


  allocate(values(npoints))
  allocate(lats(npoints))
  allocate(lons(npoints))

  call grib_open_file(ifile, modelfilename, open_mode, istatus)

  call grib_index_create(iindex, modelfilename, 'shortName,level')

  call grib_index_select(iindex, 'shortName', varname)
  call grib_index_select(iindex, 'level', level)
  call grib_new_from_index(iindex, igrib, istatus)

  call grib_get_data(igrib, lats, lons, values, istatus)

  ! test
  call mymax(values, npoints, mymax_test)
  call mymin(values, npoints, mymin_test)
  write(*,*) "shortName=", varname, "level=", level, "max=", mymax_test, "min=", mymin_test

  call grib_release(igrib, istatus)

  call grib_index_release(iindex, istatus)
  call grib_close_file(ifile, istatus)

  ! deallocate
  deallocate(lats)
  deallocate(lons)

end subroutine read_surface_field

!----------------------------------------------------------------------------

subroutine read_surface_field_latlon(modelfilename, varname, level, values, lats, lons)
  use parkind1, only : jpim, jprb
  use rttovscatt_mod, only : npoints
  use grib_api

  implicit none

#include "test_MinMax.interface"

  ! dummy arguments
  character(LEN=256), intent(in)                             :: modelfilename
  character(LEN=32), intent(in)                              :: varname
  integer(KIND=jpim), intent(in)                             :: level
  real(KIND=jprb), intent(out), dimension(:), allocatable    :: values ! (npoints)
  real(KIND=jprb), intent(out), dimension(:), allocatable    :: lats ! (npoints)
  real(KIND=jprb), intent(out), dimension(:), allocatable    :: lons ! (npoints)

  ! local variables
  character(LEN=*), parameter        :: open_mode="r"
  integer (KIND=jpim)                :: iindex
  integer (KIND=jpim)                :: ifile
  integer (KIND=jpim)                :: istatus
  integer (KIND=jpim)                :: igrib
  real (KIND=jprb)                   :: mymax_test
  real (KIND=jprb)                   :: mymin_test


  allocate(values(npoints))
  allocate(lats(npoints))
  allocate(lons(npoints))

  call grib_open_file(ifile, modelfilename, open_mode, istatus)

  call grib_index_create(iindex, modelfilename, 'shortName,level')

  call grib_index_select(iindex, 'shortName', varname)
  call grib_index_select(iindex, 'level', level)
  call grib_new_from_index(iindex, igrib, istatus)

  call grib_get_data(igrib, lats, lons, values, istatus)

  ! test
  call mymax(values, npoints, mymax_test)
  call mymin(values, npoints, mymin_test)
  write(*,*) "shortName=", varname, "level=", level, "max=", mymax_test, "min=", mymin_test

  call grib_release(igrib, istatus)

  call grib_index_release(iindex, istatus)
  call grib_close_file(ifile, istatus)

end subroutine read_surface_field_latlon
