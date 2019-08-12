subroutine read_sva_fields(ioin, svafilename, values, nprof) ! values(1,:) zenith values(2,:) azimuth
	
	use parkind1
	implicit none

	! dummy arguments
	integer(KIND=jpim), intent(in)   :: ioin
	integer(KIND=jpim), intent(in)   :: nprof
	character(LEN=256), intent(in)   :: svafilename
	real(KIND=jprb), intent(out), dimension(:, :), allocatable :: values

	! local variables
	integer(KIND=jpim)                :: ios
	integer(KIND=jpim)                :: iprof
	character(LEN=256)                :: skipstring


	allocate(values(2, nprof))

	open(ioin, file=trim(svafilename), status='old', iostat=ios)
	if (ios /= 0) then
    	WRITE(*,*) 'error opening sva file ios= ', ios
 	endif

! skip the first two lines

 	read(ioin, *) skipstring
 	read(ioin, *) skipstring

! read the data
 	DO iprof = 1, nprof
 		read(ioin, *) values(:, iprof)
 	ENDDO

 	CLOSE(ioin)

end subroutine read_sva_fields

subroutine read_unstructured_sva_fields(ioin, svafilename, values, select_table, nvalidprof) 
	! values(1,:) zenith values(2,:) azimuth
	! select_table = (2, nvalidprof) 
	
	use parkind1
	use rttovscatt_mod, only: Ni, Nj, npoints, model_res, base_nw_lat, base_nw_lon
	implicit none

	! dummy arguments
	integer(KIND=jpim), intent(in)   	:: ioin
	integer(KIND=jpim), intent(in)   	:: nvalidprof
	character(LEN=256), intent(in)   	:: svafilename
	real(KIND=jprb), intent(out), dimension(:, :), allocatable :: values
	integer(KIND=jpim), intent(out), dimension(:, :), allocatable :: select_table

	! local variables
	integer(KIND=jpim)                :: ios
	integer(KIND=jpim)                :: ivalidprof
	real(KIND=jprb)					  :: lon, lat
	integer(KIND=jpim)				  :: lon_index, lat_index	

	! open the file
	open(ioin, file=trim(svafilename), status='old', iostat=ios)

	if (ios /= 0) then
    	WRITE(*,*) 'error opening sva file ios= ', ios
 	endif

 	! WRITE(*,*) '[RTTOVSCATT_SVA]: Start reading sva!'

	allocate(values(2, nvalidprof))
	allocate(select_table(2, nvalidprof))  

! read the data ！add 0.5 to avoid truncation

 	DO ivalidprof = 1, nvalidprof
 		read(ioin, *) lon, lat
 		select_table(1, ivalidprof) = INT((lon - base_nw_lon)/model_res + 1 + 0.5)
 		select_table(2, ivalidprof) = INT((base_nw_lat - lat)/model_res + 1 + 0.5)
 		read(ioin, *) values(:, ivalidprof)
 	ENDDO
 	
 	CLOSE(ioin)

 	! WRITE(*,*) '[RTTOVSCATT_SVA]: Finish reading sva!'

end subroutine read_unstructured_sva_fields