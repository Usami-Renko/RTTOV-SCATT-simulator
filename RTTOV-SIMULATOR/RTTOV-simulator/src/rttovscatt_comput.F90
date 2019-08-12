!====================Computational routines==============================

! p(1:nlevels), ph(1:nlevels) (for one time) ph(0) = 0.
! no api, we apply it manually

subroutine generate_std_levels(p, ph)
  use parkind1, only : jprb, jpim
  implicit none

  ! dummy arguments
  real(KIND=jprb), intent(in), dimension(:)                       :: p
  real(KIND=jprb), intent(out), dimension(:), allocatable         :: ph

  ! local variables
  integer(KIND=jpim)  :: ilev
  integer(KIND=jpim)  :: nlevels

  nlevels = size(p)
  allocate(ph(nlevels))

  ph(1) = 0.
  DO ilev=2, nlevels
    !geometric average
    ph(ilev) = sqrt(p(ilev-1)*p(ilev))

    !ln average 
    !ph(i) = exp((log(p(i-1))+log(p(i)))/2)
  ENDDO

end subroutine generate_std_levels