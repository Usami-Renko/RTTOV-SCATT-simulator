PROGRAM test

    implicit none

    real :: q_ext, q_sct, q_asm, q_bsct ! [cm^2]
    real :: f,t,dmax,cabs,csca,cbsc,g,ph(37),re
    integer :: nshp,iret,is_loaded

    f    = 89.
    t    = 250. 
    nshp = 9
    dmax = 5000.
    is_loaded = 0 

    call scatdb(f,t,nshp,dmax,cabs,csca,cbsc,g,ph,re,iret,is_loaded)

    if (iret /= 0) then
        write(*,*) 'Aborting due to problem calling Liu(2008) DDA parameterisation: ', iret
        stop
    endif

    ! Convert output [m^2] cross sections to [cm^2] and return extinction
    q_ext = (csca + cabs) * 1e4
    q_sct = csca * 1e4
    q_bsct = cbsc * 1e4

    ! Output asymmetry parameter is OK (apart from JPRM kind)
    q_asm = g

    write(*,*) q_ext, q_sct, q_asm, q_bsct

END PROGRAM test