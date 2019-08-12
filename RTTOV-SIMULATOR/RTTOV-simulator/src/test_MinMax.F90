subroutine mymax(values, length, max)
      use parkind1

      implicit none
      real(KIND=jprb), dimension(:), intent(in)                  :: values
      integer (KIND=jpim), intent(in)                            :: length
      real(KIND=jprb), intent(out)                               :: max

      integer(KIND=jpim)     :: i

      max = values(1)
      do i=2, length
        if (values(i)>max) then
          max = values(i)
        end if
      end do 

end subroutine mymax

subroutine mymin(values, length, min)
      use parkind1
      
      implicit none
      real(KIND=jprb), dimension(:), intent(in)                  :: values
      integer (KIND=jpim), intent(in)                            :: length
      real(KIND=jprb), intent(out)                               :: min

      integer(KIND=jpim)     :: i

      min = values(1)
      do i=2, length
        if (values(i)<min) then
          min = values(i)
        end if
      end do 

end subroutine mymin

subroutine mymin_index(values, length, min_index)
      use parkind1
      
      implicit none
      real(KIND=jprb), dimension(:), intent(in)                  :: values
      integer (KIND=jpim), intent(in)                            :: length
      integer(KIND=jpim), intent(out)                            :: min_index

      integer(KIND=jpim)     :: i
      real(KIND=jprb)        :: min

      min = values(1)
      min_index = 1
      do i=2, length
        if (values(i)<min) then
          min = values(i)
          min_index = i
        end if
      end do 

end subroutine mymin_index