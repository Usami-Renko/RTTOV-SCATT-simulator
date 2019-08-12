
subroutine print_1d_data(oned_data, io)
	USE parkind1, only: jpim, jprb
	USE rttovscatt_mod, only: FMT1, FMT2, ncolumn 

	implicit none

	! dummy variables
	real(KIND=jprb), dimension(:), intent(in)  :: oned_data
	integer(KIND=jpim), intent(in) :: io

	! my variables
	INTEGER(KIND=jpim) :: idata, ndata, joff
  	INTEGER(KIND=jpim) :: iline, nline
  	INTEGER(KIND=jpim) :: nmod, imod

  	ndata = size(oned_data)

  	nline = ndata / ncolumn
  	nmod = ndata - nline * ncolumn
  	IF (nmod .ne. 0) THEN
  			nline = nline + 1
  	ENDIF

  	DO iline = 1, nline
  		joff = (iline-1) * ncolumn
  		IF (iline .ne. nline) THEN
  			write(io, FMT1) oned_data(joff+1:joff+ncolumn)
  		ELSE
  			DO imod = 1, nmod
  				write(io, FMT2, advance="no") oned_data(joff+imod)
  			ENDDO
  			write(io,'(a)',advance="yes")
  		ENDIF
  	ENDDO 

end subroutine print_1d_data