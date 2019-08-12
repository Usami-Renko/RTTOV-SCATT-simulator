module rttovscatt_mod
	
	use parkind1

	implicit none


	real(KIND=jprb)         		:: pi = 3.141592654

	! real(KIND=jprb)				  :: model_res = 0.03
    ! integer(KIND=jpim)              :: npoints = 4179171
	! integer(KIND=jpim)              :: Ni = 2501
	! integer(KIND=jpim)              :: Nj = 1671

	! real(KIND=jprb)                 :: base_nw_lat = 60.1
	! real(KIND=jprb)                 :: base_nw_lon = 70.0
	! real(KIND=jprb)                 :: base_se_lat = 10.0
	! real(KIND=jprb)                 :: base_se_lon = 145.0

	! real(KIND=jprb)				  :: model_res = 0.03
    ! integer(KIND=jpim)              :: npoints = 1212201
	! integer(KIND=jpim)              :: Ni = 1101
	! integer(KIND=jpim)              :: Nj = 1101

	! integer(KIND=jpim)              :: base_nw_lat = 50
	! integer(KIND=jpim)              :: base_nw_lon = 102
	! integer(KIND=jpim)              :: base_se_lat = 17
	! integer(KIND=jpim)              :: base_se_lon = 135
	
    real(KIND=jprb)				    :: model_res = 0.1
	integer(KIND=jpim)              :: npoints = 376251
	integer(KIND=jpim)              :: Ni = 751
	integer(KIND=jpim)              :: Nj = 501

	integer(KIND=jpim)              :: base_nw_lat = 65
	integer(KIND=jpim)              :: base_nw_lon = 70
	integer(KIND=jpim)              :: base_se_lat = 15
	integer(KIND=jpim)              :: base_se_lon = 145

	
	! integer(KIND=jpim)              :: level_index_start = 5 ! 1000 --> 10
	! integer(KIND=jpim)              :: level_index_stop  = 23

	! real(KIND=jprb)                 :: plevels(19) = (/ &
	! 100., 200., 300., 350., 400., 450., 500., 550., 600., 650.,&
	! 700., 750., 800., 850., 900., 925., 950., 975.,1000./)

	integer(KIND=jpim)              :: level_index_start = 4 ! 1000 --> 10
	integer(KIND=jpim)              :: level_index_stop  = 33

	real(KIND=jprb)                 :: plevels(30) = (/ &
	   10.,  20.,  30.,  50.,  70., 100., 125., 150., 175., 200.,&
	  225., 250., 275., 300., 350., 400., 450., 500., 550., 600.,&
	  650., 700., 750., 800., 850., 900., 925., 950., 975.,1000./)

	
	
	real(KIND=jprb)                 :: delta = 0.0001
	CHARACTER(LEN=32)               :: FMT1 = '(1X, 5e18.8)'
	CHARACTER(LEN=32)               :: FMT2 = '(1X, e18.8)'
	integer(KIND=jpim)              :: ncolumn = 5
	real(KIND=jprb)                 :: FillValue = -9999999.9
	

end module rttovscatt_mod
