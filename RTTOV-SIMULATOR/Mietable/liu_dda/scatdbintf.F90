SUBROUTINE SCATDBINTF(Ts, Fs, Ds, NSHP, OP, nT, nF, nD)
!
!       INPUT  : Ts, Fs, Ds: 
!                Temperatures [K], 
!                Frequencies [GHz], 
!                Dmaxs [mm]  
!                nT, nF, nD: Number of T, F, D (hide)     
!                NSHP: shape ID
!       OUTPUT : OP: optical properties: 
!                Cext [m^2], Csca [m^2], g [-]
    IMPLICIT NONE
    
    ! INTERFACE VARS
    INTEGER :: nT, nF, nD
    INTEGER :: NSHP
    REAL    :: Ts(nT), Fs(nF), Ds(nD)
    REAL    :: OP(nT, nF, nD, 3)
    ! LOCAL VARS
    INTEGER :: IR, IL
    INTEGER :: Ti, Fi, Di 
    REAL    :: T, F, D
    REAL    :: CABS, CSCA, CBSC, G, PH(37), RE

!f2py   INTEGER INTENT(HIDE), DEPEND(Ts) :: nT=SHAPE(Ts,0)
!f2py   INTEGER INTENT(HIDE), DEPEND(Fs) :: nF=SHAPE(Fs,0)
!f2py   INTEGER INTENT(HIDE), DEPEND(Ds) :: nD=SHAPE(Ds,0) 
!f2py   INTENT(OUT, HIDE) OP

    IL = 0
    IR = 0
    OP = -9999.
    Ds = Ds * 1e3 ! [mm] --> [um] 

    DO Fi = 1, nF
    DO Ti = 1, nT
    DO Di = 1, nD
        F = Fs(Fi)
        T = Ts(Ti)
        D = Ds(Di)

        ! write(*,*) F, T, D, NSHP

        CALL SCATDB(F,T,NSHP,D,CABS,CSCA,CBSC,G,PH,RE,IR,IL)
        
        IF (IR /= 0) THEN
            WRITE(*,*) 'Error ', IR
            STOP
        ENDIF
        OP(Ti, Fi, Di, 1) = CABS + CSCA
        OP(Ti, Fi, Di, 2) = CSCA
        OP(Ti, Fi, Di, 3) = G
    ENDDO
    ENDDO
    ENDDO

END SUBROUTINE SCATDBINTF
