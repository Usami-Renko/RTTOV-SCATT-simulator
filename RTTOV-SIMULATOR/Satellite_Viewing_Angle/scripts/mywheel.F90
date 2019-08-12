C FILE: MYWHEEL.F
        SUBROUTINE CHECK_TABLE(IP, IT, NLA, NLO, GT, NP)
C
C       INPUT  : IP: (NP, 2) LON/LAT IT: TYPICALY 1
C       OUTPUT : GT: INTIALIZED AS (NLA, NLO, 4) 0 CHECK AS 1     
C
        IMPLICIT NONE
        INTEGER NP, IT, NLA, NLO, ILAT, ILON, IPOINT
        INTEGER :: IP(NP, 2), GT(NLA, NLO, 4)
        INTEGER :: CHECK_EXT(4) ! DOWN UP LEFT RIGHT
        INTEGER :: NOWLAT, NOWLON

Cf2py   INTENT(IN, COPY) IP
Cf2py   INTENT(IN, COPY) IT, NLA, NLO
Cf2py   INTEGER INTENT(HIDE), DEPEND(IP) :: NP=SHAPE(IP,0) 
Cf2py   INTENT(OUT, HIDE) GT

        ! INITIALZE
        GT = 0

        DO IPOINT = 1, NP
            NOWLAT = IP(IPOINT, 2)
            NOWLON = IP(IPOINT, 1)
            
            CHECK_EXT = IT

            IF ( NOWLAT < IT ) THEN
                CHECK_EXT(1) = NOWLAT
            ENDIF

            IF ( NLA-NOWLAT < IT ) THEN
                CHECK_EXT(2) = NLA-NOWLAT
            ENDIF

            IF ( NOWLON < IT ) THEN
                CHECK_EXT(3) = NOWLON
            ENDIF

            IF ( NLO-NOWLON < IT ) THEN
                CHECK_EXT(4) = NLO-NOWLON
            ENDIF 

            ! get the four direction of GT

            DO ILAT = NOWLAT-CHECK_EXT(1)+1, NOWLAT
                DO ILON = NOWLON-CHECK_EXT(3)+1, NOWLON
                    GT(ILAT, ILON, 1) = GT(ILAT, ILON, 1)+1
                ENDDO
            ENDDO

            IF ( CHECK_EXT(2) > 0) THEN

            DO ILAT = NOWLAT+1, NOWLAT+CHECK_EXT(2)
                DO ILON = NOWLON-CHECK_EXT(3)+1, NOWLON
                    GT(ILAT, ILON, 2) = GT(ILAT, ILON, 2)+1
                ENDDO
            ENDDO

            ENDIF

            IF ( CHECK_EXT(4) > 0) THEN

            DO ILAT = NOWLAT-CHECK_EXT(1)+1, NOWLAT
                DO ILON = NOWLON+1, NOWLON+CHECK_EXT(4)
                    GT(ILAT, ILON, 3) = GT(ILAT, ILON, 3)+1
                ENDDO
            ENDDO

            ENDIF

            IF ( CHECK_EXT(4) > 0 .AND. CHECK_EXT(2) > 0 ) THEN

            DO ILAT = NOWLAT+1, NOWLAT+CHECK_EXT(2)
                DO ILON = NOWLON+1, NOWLON+CHECK_EXT(4)
                    GT(ILAT, ILON, 4) = GT(ILAT, ILON, 4)+1
                ENDDO
            ENDDO

            ENDIF

        ENDDO

        END

C END OF FILE MYWHEEL.F