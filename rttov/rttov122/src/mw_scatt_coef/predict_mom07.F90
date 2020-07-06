 subroutine predict_mom07(m2,tc,n,m)

! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.

 ! subroutine to predict nth moment m given m2, tc, n
 ! if m2=-9999 then the routine will predict m2 given m,n,tc

 ! History:
 ! Version   Date        Comment
 ! -------   ----        -------
 !   1.0     15/01/09    Copied from 
 !                      http://www.mmm.ucar.edu/people/field/files/web/Refereed%20publications.html
 !                       (Amy Doherty)
 !   1.1     27/04/10    Given an F90 interface (Alan Geer)

  use parkind1, only: jprb

  implicit none

  real (kind=jprb), intent(inout) :: m2, m
  real (kind=jprb), intent(in)    :: tc, n

!INTF_END

  real (kind=jprb) :: a1,a2,a3,b1,b2,b3,c1,c2,c3
  real (kind=jprb) :: a_,b_,c_,A,B,C

  a1=      13.6078_jprb
  a2=     -7.76062_jprb
  a3=     0.478694_jprb
  b1=   -0.0360722_jprb
  b2=    0.0150830_jprb
  b3=   0.00149453_jprb
  c1=     0.806856_jprb
  c2=   0.00581022_jprb
  c3=    0.0456723_jprb

  a_=a1+a2*n+a3*n**2.0_jprb
  b_=b1+b2*n+b3*n**2.0_jprb
  c_=c1+c2*n+c3*n**2.0_jprb

  A=exp(a_)
  B=b_
  C=c_

 ! predict m from m2 and tc
  if(m2.ne.-9999._jprb) then
    m=A*exp(B*tc)*m2**C
  else
 ! get m2 if mass-dimension relationship not proportional to D**2
    m2=(m/(A*exp(B*tc)))**(1.0_jprb/C)
  endif

  return
  end

