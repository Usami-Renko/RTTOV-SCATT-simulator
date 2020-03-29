PROGRAM test

    implicit none

    real :: q_ext, q_sct, q_asm, q_bsct ! [cm^2]
    real :: f,t,dmax,cabs,csca,cbsc,g,ph(37),re
    integer :: nshp,iret,is_loaded

    f    = 10.65
    t    = 250. 
    nshp = 9
    dmax = 10000.
    is_loaded = 0 

    call scatdb(f,t,nshp,dmax,cabs,csca,cbsc,g,ph,re,iret,is_loaded)

    if (iret /= 0) then
        write(*,*) 'Aborting due to problem calling Liu(2008) DDA parameterisation: ', iret
        stop
    endif

    q_ext = (csca + cabs)
    q_sct = csca
    q_bsct = cbsc

    ! Output asymmetry parameter is OK (apart from JPRM kind)
    q_asm = g

    write(*,*) q_ext, q_sct, q_asm, q_bsct

END PROGRAM test