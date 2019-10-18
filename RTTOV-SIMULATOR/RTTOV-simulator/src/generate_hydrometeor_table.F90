Program generate_hydrometeor_table
	USE parkind1, only: jpim, jprb

	implicit none
#include "rttovscatt_gribapi.interface"
#include "rttovscatt_select_level.interface"
#include "rttovscatt_select_levels.interface"
#include "print_1d_data.interface"
#include "test_MinMax.interface"


! io units

	INTEGER(KIND=jpim), PARAMETER :: io  = 20

 	CHARACTER(LEN=256) :: model_filename
 	CHARACTER(LEN=256) :: output_filename

! read variables

 	real(KIND=jprb), allocatable    :: levels_raw(:, :)         	!(nlevels, npoints)
  real(KIND=jprb), allocatable    :: levels_selected_raw(:, :)    !(nlevels, nprof)
  real(KIND=jprb), allocatable    :: level_raw(:)                 !(npoints)
  real(KIND=jprb), allocatable    :: level_selected_raw(:)        !(nprof)

  integer(KIND=jpim)              :: nw_lat
  integer(KIND=jpim)              :: nw_lon
  integer(KIND=jpim)              :: se_lat
  integer(KIND=jpim)              :: se_lon

  real(KIND=jprb), dimension(:), allocatable    :: lats_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lons_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lats_selected_raw ! (npoints)
  real(KIND=jprb), dimension(:), allocatable    :: lons_selected_raw ! (npoints)

  INTEGER            							  :: ios
  CHARACTER(LEN=32)                             :: shortName

  real(KIND=jprb) :: centre_lat, centre_lon
  integer(KIND=jpim) :: centre_index

! loop variables
	INTEGER(KIND=jpim) :: ilev, nlevels
	INTEGER(KIND=jpim) :: nprof

 !- End of header------------------------------------------


  !=====================================================
  !========== Interactive inputs == start ==============

 	WRITE(0,*) 'enter path of MODEL file'
  READ(*,*) model_filename
  ! naming output_filename//"_rain/_sp/_clw/_ciw/_cc/_gh"
  WRITE(0,*) 'enter path of OUTPUT file'
  READ(*,*) output_filename
  WRITE(0,*) 'enter number of profiles'
  READ(*,*) nprof
  WRITE(0,*) 'enter number of profile levels'
  READ(*,*) nlevels
  WRITE(0,*) 'enter nw lat'
  READ(*,*) nw_lat
  WRITE(0,*) 'enter nw lon'
  READ(*,*) nw_lon
  WRITE(0,*) 'enter se lat'
  READ(*,*) se_lat
  WRITE(0,*) 'enter se lon'
  READ(*,*) se_lon

  !========== Interactive inputs == end ==============
  !===================================================

  write(*, *) 'start reading profile!'

  ! generate the hydrometeor table

  shortName = "ccl"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  OPEN(io, file=TRIM(output_filename)//"_cc.dat", status='unknown', form='formatted', iostat=ios)
  DO ilev = 1, nlevels	
  	call print_1d_data(levels_selected_raw(ilev,:), io)
  ENDDO
  CLOSE(io, iostat=ios)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "clwmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  OPEN(io, file=TRIM(output_filename)//"_clw.dat", status='unknown', form='formatted', iostat=ios)
  DO ilev = 1, nlevels	
  	call print_1d_data(levels_selected_raw(ilev,:), io)
  ENDDO
  CLOSE(io, iostat=ios)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "icmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  OPEN(io, file=TRIM(output_filename)//"_ciw.dat", status='unknown', form='formatted', iostat=ios)
  DO ilev = 1, nlevels	
  	call print_1d_data(levels_selected_raw(ilev,:), io)
  ENDDO
  CLOSE(io, iostat=ios)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "rwmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  OPEN(io, file=TRIM(output_filename)//"_rain.dat", status='unknown', form='formatted', iostat=ios)
  DO ilev = 1, nlevels	
  	call print_1d_data(levels_selected_raw(ilev,:), io)
  ENDDO
  CLOSE(io, iostat=ios)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "snmr"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  OPEN(io, file=TRIM(output_filename)//"_sp.dat", status='unknown', form='formatted', iostat=ios)
  DO ilev = 1, nlevels	
  	call print_1d_data(levels_selected_raw(ilev,:), io)
  ENDDO
  CLOSE(io, iostat=ios)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  shortName = "gh"
  call read_levels_field(model_filename, shortName, levels_raw, nlevels)
  call select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  OPEN(io, file=TRIM(output_filename)//"_gh.dat", status='unknown', form='formatted', iostat=ios)
  DO ilev = 1, nlevels	
  	call print_1d_data(levels_selected_raw(ilev,:), io)
  ENDDO
  CLOSE(io, iostat=ios)
  deallocate(levels_selected_raw)
  deallocate(levels_raw)

  ! generate the storm-centre and print lat / lon

  shortName = "sp"
  call read_surface_field_latlon(model_filename, shortName, 0, level_raw, lats_raw, lons_raw)
  level_raw = level_raw / 100.0
  call select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  call select_values(lats_raw, lats_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  call select_values(lons_raw, lons_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  
  call mymin_index(level_selected_raw, nprof, centre_index)
  centre_lat = lats_selected_raw(centre_index)
  centre_lon = lons_selected_raw(centre_index)
  write(*,'(a11, f8.1, 1x, a11, f8.1)') "centre_lat=", centre_lat, "centre_lon=", centre_lon

  OPEN(io, file=TRIM(output_filename)//"_centre.dat", status='unknown', form='formatted', iostat=ios)
  write(io, '(2f8.1)') centre_lat, centre_lon
  call print_1d_data(lats_selected_raw, io)
  call print_1d_data(lons_selected_raw, io)
  CLOSE(io, iostat=ios)

  deallocate(level_raw)
  deallocate(lats_raw)
  deallocate(lons_raw)
  deallocate(level_selected_raw)
  deallocate(lats_selected_raw)
  deallocate(lons_selected_raw)
  

End Program generate_hydrometeor_table