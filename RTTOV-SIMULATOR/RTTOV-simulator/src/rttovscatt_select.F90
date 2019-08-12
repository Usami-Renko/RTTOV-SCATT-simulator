!====================select_values==============================

subroutine select_levels_values(levels_raw, levels_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  use parkind1, only : jprb, jpim
  use rttovscatt_mod, only: npoints, base_nw_lat, base_nw_lon, model_res

  implicit none

#include "rttovscatt_select_level.interface"

  !dummy arguments
  real(KIND=jprb), intent(in), dimension(:, :)                        :: levels_raw
  integer(KIND=jpim), intent(in)                                      :: nw_lat
  integer(KIND=jpim), intent(in)                                      :: nw_lon
  integer(KIND=jpim), intent(in)                                      :: se_lat
  integer(KIND=jpim), intent(in)                                      :: se_lon 
  real(KIND=jprb), intent(out), allocatable, dimension(:, :)          :: levels_selected_raw

  ! levels_raw (nlevels, npoints)
  ! local variables

  integer(KIND=jpim) :: index_lat_n
  integer(KIND=jpim) :: index_lat_s
  integer(KIND=jpim) :: index_lon_w
  integer(KIND=jpim) :: index_lon_e
  integer(KIND=jpim) :: selected_length_lon
  integer(KIND=jpim) :: selected_length_lat
  integer(KIND=jpim) :: selected_npoints
  integer(KIND=jpim) :: ilev
  integer(KIND=jpim) :: nlevels
  real(KIND=jprb), allocatable    :: level_selected_raw_tmp(:)
  real(KIND=jprb), allocatable    :: level_raw_tmp(:)  


  index_lat_n = INT((base_nw_lat - nw_lat)/model_res + 1 + 1e-6)
  index_lat_s = INT((base_nw_lat - se_lat)/model_res + 1 + 1e-6)
  index_lon_w = INT((nw_lon - base_nw_lon)/model_res + 1 + 1e-6)
  index_lon_e = INT((se_lon - base_nw_lon)/model_res + 1 + 1e-6)
  selected_length_lat = index_lat_s - index_lat_n + 1
  selected_length_lon = index_lon_e - index_lon_w + 1
  selected_npoints = selected_length_lat * selected_length_lon
  nlevels = size(levels_raw, dim=1)

  allocate(levels_selected_raw(nlevels, selected_npoints))

  allocate(level_raw_tmp(npoints))
  DO ilev = 1, nlevels   
      level_raw_tmp = levels_raw(ilev, :)
      call select_values(level_raw_tmp, level_selected_raw_tmp, nw_lat, nw_lon, se_lat, se_lon)
      levels_selected_raw(ilev, :) = level_selected_raw_tmp
      deallocate(level_selected_raw_tmp)
  ENDDO
  deallocate(level_raw_tmp)

end subroutine select_levels_values


!-------------------------------------------------------------------------------

subroutine select_values(level_raw, level_selected_raw, nw_lat, nw_lon, se_lat, se_lon)
  use parkind1, only : jprb, jpim
  use rttovscatt_mod, only: npoints, base_nw_lat, base_nw_lon, Ni, Nj, model_res

  implicit none

  !dummy arguments
  real(KIND=jprb), intent(in), dimension(:)                           :: level_raw
  integer(KIND=jpim), intent(in)                                      :: nw_lat
  integer(KIND=jpim), intent(in)                                      :: nw_lon
  integer(KIND=jpim), intent(in)                                      :: se_lat
  integer(KIND=jpim), intent(in)                                      :: se_lon 
  real(KIND=jprb), intent(out), allocatable, dimension(:)             :: level_selected_raw


  !local variables
  real(KIND=jprb), allocatable      :: level_field(:, :)
  real(KIND=jprb), allocatable      :: level_selected_field(:, :)
  integer(KIND=jpim) :: index_lat_n
  integer(KIND=jpim) :: index_lat_s
  integer(KIND=jpim) :: index_lon_w
  integer(KIND=jpim) :: index_lon_e
  integer(KIND=jpim) :: selected_length_lon
  integer(KIND=jpim) :: selected_length_lat
  integer(KIND=jpim) :: selected_npoints

  index_lat_n = INT((base_nw_lat - nw_lat)/model_res + 1 + 1e-6)
  index_lat_s = INT((base_nw_lat - se_lat)/model_res + 1 + 1e-6)
  index_lon_w = INT((nw_lon - base_nw_lon)/model_res + 1 + 1e-6)
  index_lon_e = INT((se_lon - base_nw_lon)/model_res + 1 + 1e-6)
  selected_length_lat = index_lat_s - index_lat_n + 1
  selected_length_lon = index_lon_e - index_lon_w + 1  
  selected_npoints = selected_length_lat * selected_length_lon

  allocate(level_field(Ni, Nj))       ! (fast)lon-->lat(slow)
  allocate(level_selected_field(selected_length_lon, selected_length_lat))
  allocate(level_selected_raw(selected_npoints))
  

  level_field = reshape(level_raw, (/Ni, Nj/))
  level_selected_field = level_field(index_lon_w:index_lon_e, index_lat_n:index_lat_s)
  level_selected_raw = reshape(level_selected_field, (/selected_npoints/))


  deallocate(level_field)
  deallocate(level_selected_field)

end subroutine select_values

subroutine select_levels_unstructured_values(levels_raw, levels_selected_raw, select_table, nvalidprofs)
  use parkind1, only : jprb, jpim
  use rttovscatt_mod, only: npoints

  implicit none

#include "rttovscatt_select_level.interface"

  !dummy arguments
  real(KIND=jprb), intent(in), dimension(:, :)                        :: levels_raw
  real(KIND=jprb), intent(out), allocatable, dimension(:, :)          :: levels_selected_raw
  integer(KIND=jpim), intent(in), dimension(:, :)                     :: select_table 
  integer(KIND=jpim), intent(in)                                      :: nvalidprofs

  ! levels_raw (nlevels, npoints)
  ! levels_selected_raw (nlevels, nvalidprofs)

  ! local variables

  integer(KIND=jpim) :: ilev
  integer(KIND=jpim) :: nlevels
  real(KIND=jprb), allocatable    :: level_selected_raw_tmp(:)
  real(KIND=jprb), allocatable    :: level_raw_tmp(:)  


  nlevels = size(levels_raw, dim=1)

  allocate(levels_selected_raw(nlevels, nvalidprofs))

  allocate(level_raw_tmp(npoints))
  DO ilev = 1, nlevels   
      level_raw_tmp = levels_raw(ilev, :)
      call select_unstructured_values(level_raw_tmp, level_selected_raw_tmp, select_table, nvalidprofs)
      levels_selected_raw(ilev, :) = level_selected_raw_tmp
      deallocate(level_selected_raw_tmp)
  ENDDO
  deallocate(level_raw_tmp)


end subroutine select_levels_unstructured_values

subroutine select_unstructured_values(level_raw, level_selected_raw, select_table, nvalidprofs)
! level_raw (Ni*Nj)
! level_selected_raw (nvalidprofs)
! select_table (2, nvalidprof) 
  use parkind1, only : jprb, jpim
  use rttovscatt_mod, only: Ni, Nj

  implicit none

  !dummy arguments
  real(KIND=jprb), intent(in), dimension(:)                           :: level_raw
  real(KIND=jprb), intent(out), allocatable, dimension(:)             :: level_selected_raw
  integer(KIND=jpim), intent(in), dimension(:, :)                     :: select_table 
  integer(KIND=jpim), intent(in)                                      :: nvalidprofs

  !local variables
  real(KIND=jprb), allocatable      :: level_field(:, :)
  integer(KIND=jpim)                :: ivalidprof
  integer(KIND=jpim)                :: lon_index, lat_index

  allocate(level_field(Ni, Nj))       ! (fast)lon-->lat(slow)
  allocate(level_selected_raw(nvalidprofs))

  level_field = reshape(level_raw, (/Ni, Nj/))
  
  DO ivalidprof=1, nvalidprofs
    lon_index = select_table(1, ivalidprof)
    lat_index = select_table(2, ivalidprof)
    level_selected_raw(ivalidprof) = level_field(lon_index, lat_index)
  ENDDO

  deallocate(level_field)

end subroutine select_unstructured_values

