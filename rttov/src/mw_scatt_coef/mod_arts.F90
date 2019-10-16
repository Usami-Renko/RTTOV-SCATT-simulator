module mod_arts

! Copyright:
!
!    Copyright 2017, Chalmers University of Technology/EUMETSAT, All Rights Reserved.

! Any habit from the ARTS SSP database needs to be registered here.
! That is, this module needs to be updated whenever new habits are added
! or existing ones modified.
! The order of entries in all per-habit data arrays below has to be
! identical.

! History:
! Version   Date        Comment
! -------   ----        -------
!           31/10/2017  initial version (Jana Mendrok)

use parkind1,  only: jprb, jpim, jplm
use mod_mie, only: lmax ! add whatever is needed

implicit none


!* base location of the RTTOV-SSP files from the ARTS SPP database
character(len=lmax) :: arts_folder = 'artsdb/'

!* total number of registered ARTS habits
integer (kind=jpim), parameter :: n_arts_habits = 17

!* maximum (allowed) f, T, D dimensions from ARTS-SSDB data
integer (kind=jpim), parameter :: nf_max_arts = 34
integer (kind=jpim), parameter :: nT_max_arts =  5
integer (kind=jpim), parameter :: nD_max_arts = 55

!* maximum (allowed) extrapolation factors in f, T, D from ARTS-SSDB data
real (kind=jprb), parameter :: f_extpol = 0.5
real (kind=jprb), parameter :: T_extpol = 0.5
real (kind=jprb), parameter :: D_extpol = 1e-6    ! basically 0., but we need a little
                                                  ! bit due to number representation issues

!* filenames of the habits
character(len=33), dimension (n_arts_habits) :: arts_files = (/ &       ! ID
  'PlateType1.rssp                  ', &                                !  1
  'ColumnType1.rssp                 ', &                                !  2
  '6-BulletRosette.rssp             ', &                                !  3
  'Perpendicular4-BulletRosette.rssp', &                                !  4
  'Flat3-BulletRosette.rssp         ', &                                !  5
  'IconCloudIce.rssp                ', &                                !  6
  'SectorSnowflake.rssp             ', &                                !  7
  'EvansSnowAggregate.rssp          ', &                                !  8
  '8-ColumnAggregate.rssp           ', &                                !  9
  'LargePlateAggregate.rssp         ', &                                ! 10
  'LargeColumnAggregate.rssp        ', &                                ! 11
  'LargeBlockAggregate.rssp         ', &                                ! 12
  'IconSnow.rssp                    ', &                                ! 13
  'IconHail.rssp                    ', &                                ! 14
  'GemGraupel.rssp                  ', &                                ! 15
  'LiquidSphere.rssp                ', &                                ! 16
  'BiLei_10_plates.rssp             '  &                                ! 17
  /)


!* Min / max diameters of the habits (Dmax [m])
real (kind=jprb), dimension (n_arts_habits) :: d_arts_min = (/ &
  1.315073e-05, &    !  1
  1.440998e-05, &    !  2
  1.556235e-05, &    !  3
  1.795143e-05, &    !  4
  1.987506e-05, &    !  5
  1.286199e-05, &    !  6
  2.000000e-05, &    !  7
  3.200000e-05, &    !  8
  1.942867e-05, &    !  9
  1.622964e-05, &    ! 10
  2.417353e-05, &    ! 11
  1.316186e-05, &    ! 12
  1.651748e-05, &    ! 13
  1.029416e-05, &    ! 14
  1.942867e-05, &    ! 15
  1.242944e-06, &    ! 16
  2.000000e-06  &    ! 17
  /)

real (kind=jprb), dimension (n_arts_habits) :: d_arts_max = (/ &
  1.000003e-02, &    !  1
  1.000000e-02, &    !  2
  1.000005e-02, &    !  3
  1.000003e-02, &    !  4
  1.000005e-02, &    !  5
  1.000003e-02, &    !  6
  1.023801e-02, &    !  7
  1.175533e-02, &    !  8
  9.714335e-03, &    !  9
  2.285975e-02, &    ! 10
  1.998066e-02, &    ! 11
  2.187592e-02, &    ! 12
  1.999972e-02, &    ! 13
  5.349094e-03, &    ! 14
  6.596726e-03, &    ! 15
  5.000000e-02, &    ! 16
  1.000000e-02  &    ! 17
  /)


!* Alpha and beta of mass-dimension relationship (m=a*Dmax^b) (for m [kg] and Dmax [m])
real (kind=jprb), dimension (n_arts_habits) :: alpha_arts = (/ &
  7.570440e-01, &    !  1
  3.796800e-02, &    !  2
  4.927320e-01, &    !  3
  3.248390e-01, &    !  4
  2.433350e-01, &    !  5
  1.590000e+00, &    !  6
  8.222560e-04, &    !  7
  1.963050e-01, &    !  8
  6.544800e+01, &    !  9
  2.085010e-01, &    ! 10
  2.758260e-01, &    ! 11
  3.499490e-01, &    ! 12
  3.114040e-02, &    ! 13
  3.835060e+02, &    ! 14
  1.727530e+02, &    ! 15
  5.235990e+02, &    ! 16
  2.043544e+01  &    ! 17
  /)

real (kind=jprb), dimension (n_arts_habits) :: beta_arts = (/ &
  2.477030e+00, &    !  1
  2.051090e+00, &    !  2
  2.427790e+00, &    !  3
  2.425930e+00, &    !  4
  2.425730e+00, &    !  5
  2.560000e+00, &    !  6
  1.444640e+00, &    !  7
  2.386110e+00, &    !  8
  3.000000e+00, &    !  9
  2.257080e+00, &    ! 10
  2.444020e+00, &    ! 11
  2.265690e+00, &    ! 12
  1.948600e+00, &    ! 13
  2.994190e+00, &    ! 14
  2.964610e+00, &    ! 15
  3.000000e+00, &    ! 16
  3.000000e+00  &    ! 17
  /)


end module mod_arts
