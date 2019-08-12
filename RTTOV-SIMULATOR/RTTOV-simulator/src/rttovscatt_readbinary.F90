! internally generated : FASTEM, Salinity, Surface type, watertype /
!                        Elevation, latitude, longitude

subroutine read_levels_field(iomd, modelfilename, varname, b_rects, values, nlevels)
  use parkind1, only: jpim, jprb
  use rttovscatt_mod, only: Ni, Nj, npoints
  implicit none

#include "test_MinMax.interface"

  ! dummy arguments
  integer(KIND=jpim), intent(in)                             :: iomd
  character(LEN=256), intent(in)                             :: modelfilename
  character(LEN=32), intent(in)                              :: varname
  integer(KIND=jpim), intent(in)                             :: b_rects
  real(KIND=jprb), intent(out), dimension(:, :), allocatable :: values
  integer(KIND=jpim), intent(in)                             :: nlevels

  ! local variables
  real(KIND=jprb), allocatable      :: values_table(:, :, :)
  real, allocatable                 :: read_values_temp(:, :)
  real(KIND=jprb), allocatable      :: values_temp(:, :)
  integer(KIND=jpim)                :: i, j, iz 
  real(KIND=jprb)                   :: mymax_test, mymin_test

  ! allocate local variables
  allocate(values(nlevels, npoints))
  allocate(values_table(nlevels, Ni, Nj))

  ! (nlevels: (high=>low), latitude(high=>low), longitude(low=>high))
  allocate(values_temp(Ni, Nj))
  allocate(read_values_temp(Ni, Nj))

  ! open binary file
  open(iomd, FILE=modelfilename, FORM="unformatted", ACCESS="direct", RECL=npoints)

  DO iz=1, nlevels
    ! y start form low latitude to high latitude 
    READ(iomd, rec=b_rects+iz) ((read_values_temp(i, Nj-j+1), i=1,Ni), j=1, Nj)
    values_temp = read_values_temp ! type transfer
    ! test
    call mymax(reshape(values_temp, (/npoints/)), npoints, mymax_test)
    call mymin(reshape(values_temp, (/npoints/)), npoints, mymin_test)
    write(*,'("varname=", a6, 2x, "level=", i2, 2x, "max=", e12.5, 2x, "min=", e12.5)') TRIM(varname), iz, mymax_test, mymin_test
    ! assign the value  
    values_table(nlevels+1-iz, :, :) = values_temp
  ENDDO
  
  ! squeeze the table
  values = reshape(values_table, (/nlevels, npoints/)) 

  close(iomd)

  deallocate(values_temp)
  deallocate(read_values_temp)
  deallocate(values_table)

end subroutine read_levels_field

subroutine read_surface_field(iomd, modelfilename, varname, b_rects, values)
  use parkind1, only : jpim, jprb
  use rttovscatt_mod, only : npoints, Ni, Nj
  implicit none

#include "test_MinMax.interface"

  ! dummy arguments
  integer(KIND=jpim), intent(in)                             :: iomd
  character(LEN=256), intent(in)                             :: modelfilename
  character(LEN=32), intent(in)                              :: varname
  integer(KIND=jpim), intent(in)                             :: b_rects
  real(KIND=jprb), intent(out), dimension(:), allocatable    :: values ! (npoints)

  ! local variables
  real (KIND=jprb)                   :: mymax_test
  real (KIND=jprb)                   :: mymin_test
  real, allocatable                  :: read_values_table(:, :)
  real (KIND=jprb), allocatable      :: values_table(:, :) 
  integer(KIND=jpim)                 :: i, j 


  allocate(values(npoints))

  allocate(read_values_table(Ni, Nj))
  allocate(values_table(Ni, Nj))

  open(iomd, FILE=modelfilename, FORM="unformatted", ACCESS="direct", RECL=npoints)

  ! y start form low latitude to high latitude 
  READ(iomd, rec=b_rects+1) ((read_values_table(i, Nj-j+1), i=1,Ni), j=1, Nj)

  values_table = read_values_table

  ! test
  call mymax(reshape(values_table, (/npoints/)), npoints, mymax_test)
  call mymin(reshape(values_table, (/npoints/)), npoints, mymin_test)
  write(*,'("varname=", a6, 2x, "max=", e12.5, 2x, "min=", e12.5)') TRIM(varname), mymax_test, mymin_test
    
  values = reshape(values_table, (/npoints/))

  close(iomd)

  ! deallocate
  deallocate(values_table)
  deallocate(read_values_table)

end subroutine read_surface_field