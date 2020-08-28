      PROGRAM DATAINH
!-----------------------------------------------------------------------
!     Calls CALPET, MYHRLY, POLYKH, RETENT, RLD, THERMK, TIMEX, and
!     WELCOME
!
!     Reads the data file *.INP, checks the data and performs 
!     some calculations, then writes the necessary information to 
!     the binary output file *.BIN, which serves as the input 
!     file to UNSAT-H
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      FUNVD(TD) = EXP(46.440973-(6790.4985/TD)-6.02808*LOG(TD))
      INCLUDE 'init.inc'
      DATA DV/3.01/
      LHV0 = LHV0*WATDEN
      LOG10E = LOG10(EXP(1.D0))
      CALL TIMEX
!
!     IFILE = Input filename
!
      CALL WELCOME('DATAINH',DV)
20    WRITE(LUS,1010) 
      READ(LUR,'(A)') TMFILE
      IF (INDEX(TMFILE(1:1),'0') .GT. 0) GOTO 999
      NCHR = INDEX(TMFILE,' ')
      TMFILE(NCHR:NCHR+3) = '.inp'
      OPEN(UNIT=LUI,FILE=TMFILE,STATUS='OLD',IOSTAT=INERR)
      IF (INERR .GT. 0) THEN
        WRITE(LUS,6005) IFILE
        GOTO 20
      ENDIF
6005  FORMAT(/,' No record of file ',A50,/,' Try again...',/)
      INQUIRE(UNIT=LUI,NAME=IFILE)
      WRITE(LUS,1020)
      READ(LUR,'(I)') IPRINT
      IF (IPRINT .GT. 0) THEN
        WRITE(LUS,1030)
        READ(LUR,'(I)') IOUT
        IF (IOUT .EQ. 0) THEN
          LOUT = LUS
        ELSE 
          LOUT = LUW
          LUSLUW = .TRUE.
          TMFILE(NCHR:NCHR+3) = '.lis'
          OPEN (UNIT=LOUT,FILE=TMFILE,STATUS='UNKNOWN')
          CLOSE(UNIT=LOUT,STATUS='DELETE')
          OPEN(UNIT=LOUT,FILE=TMFILE,STATUS='NEW')
          INQUIRE(UNIT=LOUT,NAME=OFILE)
          WRITE(LUS,1040) OFILE
        ENDIF
        WRITE(LOUT,1050) DV,IFILE,SDATE,STIME
      ENDIF
1010  FORMAT(/,
     &' Enter input filename without the ".INP" extension',/,           &
     &' (a "0" terminates the program) ===> ',$)
1020  FORMAT(/,                                                         &
     &' Output options in addition to creating < >.bin:',/,             &
     &'   (0) = no output',/,                                           &
     &'   (1) = options, nodes, soil properties, and initial',/,        &
     &'         conditions',/,                                          &
     &'   (2) = (1) plus PET and rain information',/,                   &
     &'   (3) = (2) plus meteorological, leaf, and bottom',/,           &
     &'             boundary flux information',//,                      &
     &' Enter choice ===> ',$)
1030  FORMAT(/,                                                         &
     &' Options for output device:',/,                                  &
     &'     (0) screen',/,                                              &
     &'     (1) file',/,                                                &
     &' Choose output device option ===> ',$)
1040  FORMAT(/,' The input summary filename is:',/,1X,A79)
1050  FORMAT(80('-'),/,"!",T30,' Program DATAINH',T80,"!",/,            &
     &"!",T32,' Version ',F4.2,T80,"!",/,80('-'),/,                     &
     &' Input Filename: ',A60,/,                                        &
     &' Date Processed: ',A11,/,' Time Processed: ',A11)
!**********************************************************************!
!     OPTIONS, CONSTANTS, AND LIMITS
!**********************************************************************!
!
!     Read TITLE
!
      READ(LUI,'(A)') TITLE
      IF (IPRINT .GT. 0) THEN
        WRITE(LOUT,'(A,/,1X,A79)') ' Title:',TITLE
        WRITE(LOUT,'(/,A,/)') ' General options:'
      ENDIF
!
!     Read IPLANT,NGRAV
!
      READ(LUI,'(2I)') IPLANT,NGRAV
      IF (IPRINT .GT. 0) WRITE(LOUT,2010) IPLANT,NGRAV
2010  FORMAT(T2,'IPLANT =',I5,T22,' NGRAV =',I5)
      G = NGRAV
!
!     Read IFDEND,IDTBEG,IDTEND
!
      READ(LUI,'(3I)') IFDEND,IDTBEG,IDTEND
      IF (IPRINT .GT. 0) WRITE(LOUT,2020) IFDEND,IDTBEG,IDTEND
      IF (IDTBEG .GT. IDTEND .OR. IDTBEG .LT. 1 .OR.                    &
     &  IDTEND .GT. 366) THEN
        WRITE(*,'(A)') ' IDTBEG and IDTEND must be between 1 and 366'
        WRITE(*,'(A)') ' and IDTBEG must be < or = than IDTEND.'
        GOTO 999
      ENDIF
2020  FORMAT(T2,'IFDEND =',I5,T22,'IDTBEG =',I5,T42,'IDTEND =',I5)
!
!     Read IYS,NYEARS,ISTEAD,IFLIST,NFLIST
!
      READ(LUI,'(5I)') IYS,NYEARS,ISTEAD,IFLIST,NFLIST
      IF (IYS .LT. 1) THEN
        WRITE(*,'(A)') ' IYS cannot be < 1. Set it to 1 for single year'
        GOTO 999
      ENDIF
      ILEAP = 0
      IF (MOD(IYS,4) .EQ. 0) THEN
        ILEAP = 1
        IF (IYS .GE. 100 .AND. MOD(IYS,100) .EQ. 0) ILEAP = 0
        IF (IYS .GE. 400 .AND. MOD(IYS,400) .EQ. 0) ILEAP = 1 
      ENDIF
      IF (IFLIST .GT. 1 .AND. NYEARS .GT. NFLIST) THEN
        WRITE(*,'(A)') ' NYEARS cannot be > NFLIST.'
        GOTO 999
      ENDIF
      IF (IPRINT .GT. 0) WRITE(LOUT,2031)IYS,NYEARS,ISTEAD,IFLIST,NFLIST
2031  FORMAT(T2,'   IYS =',I5,T22,'NYEARS =',I5,T42,'ISTEAD =',I5,/,    &
     &T2,'IFLIST =',I5,T22,'NFLIST =',I5)
      MYEAR = IYS+NYEARS-1
      IF ((MYEAR/4)-INT(MYEAR/4) .NE. 0) THEN
        IF (IDTEND .GT. 365 .OR. IFDEND .GT. 365) THEN
          WRITE(*,'(A)')                                                &
     &    ' IDTEND and/or IFDEND exceeds no. of days in last year'
          GOTO 999
        ENDIF
      ENDIF
      IF ((NYEARS .EQ. 1 .AND. IFDEND .LT. IDTBEG) .OR.                 &
     &  IFDEND .GT. IDTEND) THEN
        WRITE(*,'(A)')                                                  &
     &  ' IFDEND is not within range of IDTBEG and IDTEND.'
        GOTO 999
      ENDIF
      IF (NYEARS .EQ. 1) THEN
        NDAYS = IDTEND - IDTBEG+1
        IDEND = IFDEND
      ELSE
        NDAYS = (365+ILEAP)-IDTBEG+1
        IDEND = 365+ILEAP
      ENDIF
      IF (ISTEAD .EQ. 1 .AND. IFLIST .EQ. 3) THEN
        WRITE(*,'(A)')                                                  &
     &  ' ISTEAD = 1 option and IFLIST = 3 option are not compatible'
        GOTO 999
      ENDIF
!
!     Read NPRINT,STOPHR
!
      READ(LUI,'(I,D)') NPRINT,STOPHR
      IF (NPRINT .EQ. 1 .AND. STOPHR .EQ. 0.0) THEN
        WRITE(LUS,'(A)') ' STOPHR must be > 0.0 when NPRINT = 1'
        GOTO 999
      ENDIF
      IF (IPRINT .GT. 0) WRITE(LOUT,2021) NPRINT,STOPHR
2021  FORMAT(T2,'NPRINT =',I5,T22,'STOPHR =',1PE10.3)
!
!     Read ISMETH,INMAX,ISWDIF,DMAXBA
!
      READ(LUI,'(3I,D)') ISMETH,INMAX,ISWDIF,DMAXBA
      IF (INMAX .LT. 2) INMAX = 2
      IF (IPRINT .GT. 0) WRITE(LOUT,2071) ISMETH,INMAX,ISWDIF,DMAXBA
2071  FORMAT(T2,'ISMETH =',I5,T22,'INMAX  =',I5,T42,'ISWDIF =',I5,      &
     &T62'DMAXBA =',E10.3)
!
!     Read DELMAX,DELMIN,OUTTIM
!
      READ (LUI,'(3D)') DELMAX,DELMIN,OUTTIM
      IF (OUTTIM .EQ. 0) OUTTIM = DELMAX
      IF (OUTTIM .GT. 24.) OUTTIM = 24.0
      IF (NPRINT .EQ. 1) THEN
        IF (OUTTIM .GT. 1.) THEN
          TEST1 = 24.0-OUTTIM*NINT(24.0/OUTTIM)
          IF (TEST1 .NE. 0) THEN
            WRITE(LUS,'(A)') ' OUTTIM does not divide 24 hours evenly'
            GOTO 999
          ENDIF
        ELSE
          TEST2 = 1.0-OUTTIM*NINT(1.0/OUTTIM)
          IF (TEST2 .NE. 0) THEN
            WRITE(LUS,'(A)')' WARNING: OUTTIM does not divide 1 hour'
            WRITE(LUS,'(A)')' evenly. Problems occur if using hourly'
            WRITE(LUS,'(A)')' precip. and/or potential evaporation.'
          ENDIF
        ENDIF
      ENDIF
      IF (DELMAX .GT. OUTTIM)                                           &
     &  WRITE(LUS,'(A)') ' FYI: DELMAX is greater than OUTTIM'
      IF (IPRINT .GT. 0) WRITE(LOUT,2060) DELMAX,DELMIN,OUTTIM
2060  FORMAT(T2,'DELMAX =',1PE10.3,T22,'DELMIN =',E10.3,T42,            &
     &'OUTTIM =',E10.3)
!
!     Read RFACT,RAINIF,DHTOL,DHMAX,DHFACT
!
      READ(LUI,'(5D)') RFACT,RAINIF,DHTOL,DHMAX,DHFACT
      IF (IPRINT .GT. 0)                                                &
     &  WRITE(LOUT,2070) RFACT,RAINIF,DHTOL,DHMAX,DHFACT
2070  FORMAT(T2,' RFACT =',1PE10.3,T22,'RAINIF =',E10.3,                &
     &T42,' DHTOL =',E10.3,/,T2,' DHMAX =',E10.3,T22,'DHFACT =',E10.3)
!
!     Read KOPT,KEST,WTF
!
      READ(LUI,'(2I,D)') KOPT,KEST,WTF
      IF (IPRINT .GT. 0) WRITE(LOUT,2040) KOPT,KEST,WTF
2040  FORMAT(T2,'  KOPT =',I5,T22,'  KEST =',I5,T42,'   WTF =',1PE10.3)
      UP   = WTF
      DOWN = 1.0-UP
!
!     Read ITOPBC,IEVOPT,NFHOUR,LOWER
!
      READ(LUI,'(4I)') ITOPBC,IEVOPT,NFHOUR,LOWER
      IF (IPRINT .GT. 0) WRITE(LOUT,2009) ITOPBC,IEVOPT,NFHOUR,LOWER
2009  FORMAT(T2,'ITOPBC =',I5,T22,'IEVOPT =',I5,T42,'NFHOUR =',I5,      &
     &T62,' LOWER =',I5)
!
!     Read HIRRI,HDRY,HTOP,RHA
!
      READ(LUI,'(4D)') HIRRI,HDRY,HTOP,RHA
      IF (IPRINT .GT. 0) WRITE(LOUT,2050) HIRRI,HDRY,HTOP,RHA
2050  FORMAT(T2,' HIRRI =',1PE10.3,T22,'  HDRY =',E10.3,                &
     &T42,'  HTOP =',E10.3,T62'   RHA =',E10.3)
!
!     Read IETOPT,ICLOUD,ISHOPT
!
      READ(LUI,'(3I)') IETOPT,ICLOUD,ISHOPT
      IF (IPRINT .GT. 0) WRITE(LOUT,2030) IETOPT,ICLOUD,ISHOPT
2030  FORMAT(T2,'IETOPT =',I5,T22,'ICLOUD =',I5,T42,'ISHOPT =',I5)
!
!     Read IRAIN,HPR
!
      READ(LUI,'(I,D)') IRAIN,TMPHPR
      IF (TMPHPR .GT. 0.) HPR = TMPHPR
      IF (IPRINT .GT. 0) WRITE(LOUT,2072) IRAIN,HPR
2072  FORMAT(T2,' IRAIN =',I5,T22,'   HPR =',1PE10.3)
!
!     Read IHYS,AIRTOL,HYSTOL,HYSMXH and hysteresis restart filename
!
      IF (IPRINT .GT. 0) WRITE(LOUT,'(/,A,/)') ' Hysteresis options:'
      READ(LUI,'(I,3D,A)') IHYS,AIRTOL,HYSTOL,HYSMXH,HYFILE
      IF (IPRINT .GT. 0) THEN
        WRITE(LOUT,2041) IHYS,AIRTOL,HYSTOL,HYSMXH
        IF (IHYS .EQ. 4) THEN
!          IF (HYFILE(1:1) .EQ. '') THEN
!            WRITE(LOUT,'(A)') ' Restart filename not provided'
!            GOTO 999 
!          ELSE
             WRITE(LOUT,'(2A)') ' Restart filename = ',HYFILE
!          ENDIF
        ENDIF
      ENDIF
      IF (IHYS .GT. 0) THEN
        IF (KOPT .NE. 4) THEN
          WRITE(LUS,'(A)') ' KOPT must be 4 if hysteresis is simulated'
          GOTO 999
        ENDIF
        IF (HYSMXH .GE. HDRY) THEN
          WRITE(LUS,'(A,D)') ' Simulation stopped because HYSMXH',HYSMXH
          WRITE(LUS,'(A)') ' is not less than HDRY',HDRY
          GOTO 999
        ENDIF
      ENDIF
2041  FORMAT(T2,'  IHYS =',I5,T22,'AIRTOL =',1PE10.3,                   &
     &T42,'HYSTOL =',E10.3,T62,'HYSMXH =',E10.3)
!
!     Read IHEAT,ICONVH,DMAXHE
!
      IF (IPRINT .GT. 0) WRITE(LOUT,'(/,A,/)') ' Heat flow options:'
      READ(LUI,'(2I,D)') IHEAT,ICONVH,DMAXHE
      IF (IPRINT .GT. 0) WRITE(LOUT,2011) IHEAT,ICONVH,DMAXHE
2011  FORMAT(T2' IHEAT =',I5,T22,'ICONVH =',I5,T42,'DMAXHE =',1PE10.3)
      IF (IPLANT .EQ. 1 .AND. IHEAT .EQ. 1) THEN
        WRITE(LUS,2034) 
        IF (LUSLUW) WRITE(LUW,2034)
        GOTO 999
      ENDIF
2034  FORMAT(/,' ERROR: IHEAT not allowed to be "1" while IPLANT = 1',/)
!
!     Read UPPERH,TSMEAN,TSAMP,QHCTOP
!
      READ(LUI,'(I,3D)') UPPERH,TSMEAN,TSAMP,QHCTOP
      IF (IPRINT .GT. 0) WRITE(LOUT,2067) UPPERH,TSMEAN,TSAMP,QHCTOP
2067  FORMAT(T2,'UPPERH =',I5,T22,'TSMEAN =',1PE10.3,                   &
     &T42,' TSAMP =',E10.3,T62,'QHCTOP =',E10.3)
      IF (IETOPT .EQ. 0 .AND. (IHEAT .EQ. 1 .AND. UPPERH .EQ. 0)) THEN
        WRITE(LUS,2032) 
        IF (LUSLUW) WRITE(LUW,2032)
        GOTO 999
      ENDIF
2032  FORMAT(/,' ERROR:  IETOPT not allowed to be "0" while',/,         &
     &' IHEAT = 1 and UPPERH = 0 (surface temperature calculated)',/)
!
!     Read LOWERH,QHLEAK,TGRAD
!
      READ(LUI,'(I,2D)') LOWERH,QHLEAK,TGRAD
      IF (LOWERH .NE. 1 .AND. TGRAD .NE. 0.0) THEN
        WRITE(LOUT,2080) 
        IF (LUSLUW) WRITE(LUW,2080)
        TGRAD = 0.0
      ENDIF
2080  FORMAT(/,                                                         &
     &' NOTE:  TGRAD was reset to 0.0 because LOWERH did not equal 1',/)
      IF (IPRINT .GT. 0) WRITE(LOUT,2069) LOWERH,QHLEAK,TGRAD
2069  FORMAT(T2,'LOWERH =',I5,T22,'QHLEAK =',1PE10.3,                   &
     &T42,' TGRAD =',E10.3)
!
!     Read IVAPOR,TORT,TSOIL,VAPDIF
!
!     TORT is pore tortuosity, which is usually repor
!     TSOIL is the average soil temperature for the simulation period.
!     For Hanford, the average annual soil temperature at the 90 cm
!     depth is 15.3 C (288.46 K) for the years 1952-1980.
!     VAPDIF is the diffusivity of water vapor in air,
!     approximately 0.24 cm2/s
!
      IF (IPRINT .GT. 0) WRITE(LOUT,'(/,A,/)') ' Vapor flow options:'
      READ(LUI,'(I,3D)') IVAPOR,TORT,TSOIL,VAPDIF
      IF (IPRINT .GT. 0) WRITE(LOUT,2061) IVAPOR,TORT,TSOIL,VAPDIF
2061  FORMAT(T2,'IVAPOR =',I5,T22,'  TORT =',1PE10.3,T42,               &
     &' TSOIL =',E10.3,T62,'VAPDIF =',E10.3)
      MGR = 1.E2*MOLAR*GRAV/GASCON
      IF (IVAPOR .EQ. 1) THEN
        VAPDEN = FUNVD(TSOIL)
        VC = 3600.0*VAPDIF*TORT/WATDEN
        IF (IPRINT .GT. 3)                                              &
     &      WRITE(LOUT,4340) VAPDEN,MGR,VC
      ENDIF
4340  FORMAT(/,' IVAPOR = 1:  This option allows vapor flow',1P,/,      &
     &        ' Saturated vapor density (g/cm3)',/,                     &
     &        '  of soil when soil temperature is',/,                   &
     &        '  a constant equal to TSOIL        = ',E10.3,/,          &
     &        ' 100*MOLAR*GRAV/GASCON (K/cm)      = ',E10.3,/,          &
     &        ' VC (cm5/g/h)                      = ',E10.3)
!
!     When vapor flow is isothermal, VC is modified
!
      IF (IVAPOR .EQ. 1 .AND. IHEAT .EQ. 0) THEN
        SVD = FUNVD(TSOIL)
        VC  = VC*SVD*MGR/TSOIL
      ENDIF
!-----------------------------------------------------------------------
!     NODE INFORMATION:  Define geometric details of the problem
!-----------------------------------------------------------------------
!
!     Read MATN,NPT
!
      IF (IPRINT .GT. 0) WRITE(LOUT,'(/,A,/)') ' Grid options:'
      READ(LUI,'(2I)') MATN,NPT
      IF (IPRINT .GT. 0) WRITE(LOUT,3000) MATN,NPT
3000  FORMAT(T2,'  MATN =',I5,T22,'   NPT =',I5)
      IF (NPT .GT. M1) THEN
        WRITE(LUS,'(A)') ' NPT exceeds node dimensions M1'
        GOTO 999
      ENDIF
      IF (MATN .GT. M2) THEN
        WRITE(LUS,'(A)') ' MATN exceeds soil dimension M2' 
        GOTO 999
      ENDIF
!
!     Read MAT,Z:  Material type and depth below the surface
!
      READ(LUI,'(4(I,D))') (MAT(I), Z(I),I=1,NPT)
      DO I=1,NPT-1
        ZZ(I) = Z(I+1)-Z(I)
        IF (ZZ(I) .LE. 0.0) THEN
          WRITE(LUS,'(A,I)') ' Z does not increase at node ',I+1
          GOTO 999
        ENDIF
      ENDDO
      ZZ(NPT) = 0.0
!-----------------------------------------------------------------------
!     SOIL PROPERTY DESCRIPTION
!-----------------------------------------------------------------------
!     Read soil hydraulic property parameters
!
      CALL SHPPAR(IPRINT,LUSLUW,MATN)
!
!     Calculate water content, conductivity, and capacity at the
!     upper and lower head limits 
!
      H(1) = HIRRI
      H(2) = HDRY
      MAT2 = MAT(2)
      MAT(2) = MAT(1)
      CALL RETENT(1,2,H,1,THETA,1,C)
      CALL POLYKH(1,2,H,KL)
      SATTH = THETA(1)
      DRYTH = THETA(2)
      SATK  = KL(1)
      DRYK  = KL(2)
      SATC  = C(1)
      DRYC  = C(2)
      IF (IPRINT .GT. 0)                                                &
     &  WRITE(LOUT,4300) HIRRI,SATTH,SATK,SATC,HDRY,DRYTH,DRYK,DRYC
4300  FORMAT(80("-"),/,' Surface node bounding values:',//,             &
     &3X,'HIRRI =',1PE10.3,3X,'THETA =',E10.3,3X,'K =',E10.3,           &
     &3X,'C =',E10.3,/,4X'HDRY =',E10.3,3X,'THETA =',                   &
     &E10.3,3X,'K =',E10.3,3X,'C =',E10.3)
      MAT(2) = MAT2
!
!     Read the soil parameters for hysteresis
!
      IF (IHYS .GT. 0) CALL HYSPAR(IPRINT,LUSLUW,MATN)
!
!     Read the soil parameters for heat flow
!
      IF (IHEAT .GT. 0) CALL HEATPAR(IPRINT,LUSLUW,MATN)
!-----------------------------------------------------------------------
!     INITIAL CONDITIONS
!-----------------------------------------------------------------------
!     Input initial conditions
!     Output node, material, depth, and initial condition status
!
      IF (IPRINT .GT. 0) WRITE(LOUT,'(80("-"),/,A)')                    &
     &' Initial Conditions:'
      READ(LUI,'(I)') NDAY
      READ(LUI,'(4D)') (H(I),I=1,NPT)
      IF (IHYS .EQ. 0) THEN
        CALL RETENT(1,NPT,H,1,THETA,1,C)
        CALL POLYKH(1,NPT,H,KL)
      ELSE
        DO I=1,NPT
          HH(I) = H(I)
        ENDDO
        CALL HYSINI
        CALL HYSSHP(1,NPT,0,KL,C,THETA)
      ENDIF
      IF (IHEAT .EQ. 0) THEN
        DO I=1,NPT
          T(I) = TSOIL
        ENDDO
      ELSE
        READ(LUI,'(4D)') (T(I),I=1,NPT)
      ENDIF
      IF (IPRINT .GT. 0) THEN
        WRITE(LOUT,5010) NDAY
        WRITE(LOUT,5020) (I,Z(I),MAT(I),H(I),KL(I),C(I),THETA(I),       &
     &      T(I),I=1,NPT)
      ENDIF
5010  FORMAT(/,'   NDAY = ',I3,//,T2,'NODE',T11,'Z',T16,                &
     & 'MAT',T24,'HEAD',T35,'CONDUCTIVITY',T52,'CAPACITY',T66,'THETA',  &
     & T74,'TEMP',/,T2,'----',T9,'-----',T16,'---',T20,'------------',  &
     & T35,'------------',T50,'------------',T65,'------',T73,'-----')
5020  FORMAT(I4,F9.2,I4,1P,T20,E12.4,T35,E12.4,T50,E12.4,T65,           &
     &    0P,F6.4,F7.1)
!
      THICK(1) = 0.5*(Z(2)-Z(1))
      THICK(NPT) = 0.5*(Z(NPT)-Z(NPT-1))
      DO I=2,NPT-1
        THICK(I) = 0.5*(Z(I+1)-Z(I-1))
      ENDDO
!
      IF (IVAPOR .EQ. 1) CALL VOLVAP(1,NPT)
      TMOIST = (THETAV(1)+THETA(1))*THICK(1)                            &
     &         +(THETAV(NPT)+THETA(NPT))*THICK(NPT)
      DO I=2,NPT-1
        TMOIST = TMOIST+(THETAV(I)+THETA(I))*THICK(I)
      ENDDO
      IF (IPRINT .GT. 0) WRITE(LOUT,5030) TMOIST
5030  FORMAT(/,' Total Initial Storage (cm) = ',1PE13.6)
      IF (IHEAT .EQ. 1) THEN
        CALL THERMK(1,NPT)
        HSTORE = 0.0
        DO I=1,NPT
          HSTORE = HSTORE+CH(I)*(T(I)-T0)*THICK(I)
          IF (IVAPOR .EQ. 1) HSTORE = HSTORE+LHV0*THETAV(I)*THICK(I)
        ENDDO
!
!     convert J/cm2 to J/m2
!
        HSTORE = HSTORE*10000.
        IF (IPRINT .GT. 2) THEN
          DO I=1,NPT
            TTHETA(I) = THETA(I)
          ENDDO
          WRITE(LOUT,5024) NDAY
          WRITE(LOUT,5025) (I,Z(I),MAT(I),T(I),KH(I),CHSOIL(I),CH(I),   &
     &      I=1,NPT)
        ENDIF
        IF (IPRINT .GT. 0) WRITE(LOUT,5031) HSTORE,T0
      ENDIF
5024  FORMAT(/,' NDAY = ',I3,/,T55,'SOIL',T70,'TOTAL',/,T39,'THERMAL',  &
     &T55,'HEAT',T70,'HEAT',/,T2,'NODE',T11,'Z',T16,'MAT',T25,'TEMP',   &
     &T36,'CONDUCTIVITY',T53,'CAPACITY',T68,'CAPACITY',/,T9,'(cm)',T25, &
     &    '(K)',T37,'(J/cm/hr/K)',T53,'(J/cm3/K)',T67,'(J/cm3/K)',/,T2, &
     &    '----',T9,'-----',T16,'---',T22,'----------',T36,             &
     &    '------------',T52,'----------',T67,'----------')
5025  FORMAT(I4,F9.2,I4,1P,T20,E12.4,3E15.4,0P)
5031  FORMAT(/,' Total Initial Heat Storage = ',1PE13.7,' J/m2',        &
     &' for T0 = '3PE11.5,' K')
      IF (NPRINT .EQ. 1 .AND. IDEND-NDAY .GT. 15) WRITE(LUS,2015) IDEND
2015  FORMAT(/,' WARNING:  The choice of hourly output (NPRINT = 1)',/, &
     &' and IDEND = ',I3,' will result in a large *.RES file',/,        &
     &' when UNSAT-H is run.  Check available disk space and',/,        &
     &' reset IDEND accordingly.')
!
!     Read plant parameters
!
      IF (IPLANT .EQ. 1) CALL PLANTIN(IPRINT,LUSLUW,MATN)
!-----------------------------------------------------------------------
!     BOUNDARY CONDITIONS
!-----------------------------------------------------------------------
!
!     Read parameters for PET/MET calculations
!
      CALL ETPAR(IPRINT,LUSLUW,NFHOUR)
!
!     Read lower boundary conditions
!
      CALL LOWBC(IPRINT,LUSLUW)
!
!     Read PET/MET data
!
      IF (IFLIST .GT. 0) THEN
        CALL FILBLDR(IFLIST,NFLIST,IPRINT,ISTEAD)
        CLOSE(UNIT=LUI)
        OPEN(UNIT=LUI,FILE=PMFN(1),STATUS='OLD')
      ENDIF

      IF (IPLANT .EQ. 1 .OR. IEVOPT .GT. 0)                             &
     &  CALL ETBC(IPRINT,LUSLUW,LUI)
!
!     Partition PET values into PEVAPO and PTRANS
!
      IF (IPLANT .EQ. 1) CALL PETPART(IPRINT,LUSLUW)
!
!     Write binary file (*.BIN) for input to UNSAT-H
!
      BIFILE = TMFILE
      BIFILE(NCHR:NCHR+3) = '.bin'
      OPEN (UNIT=LUB,FILE=BIFILE,STATUS='UNKNOWN')
      CLOSE(UNIT=LUB,STATUS='DELETE')
      OPEN(UNIT=LUB,FILE=BIFILE,STATUS='NEW',FORM='UNFORMATTED')
      WRITE(LUB) DV,IPLANT,LOWER,NDAYS,NDAY,NPRINT,ITOPBC,ICONVH,MAT,   &
     &           KOPT,KEST,IVAPOR,IDEND,IFDEND,NYEARS,IYS,ISWDIF,ICLOUD,&
     &           NPT,G,MAXPOL,MAXCOE,IEVOPT,NFPET,NSOW,NHRVST,IDTEND,   &
     &           INC,MXROOT,IETOPT,ISHOPT,INMAX,ISTEAD,ILEAP,           &
     &           IRAIN,IHEAT,UPPERH,LOWERH,IFILE,HYFILE,TITLE
      WRITE(LUB) DMAXBA,DELMAX,DELMIN,RAINIF,RFACT,HIRRI,HDRY,SATK,     &
     &           DRYK,SATC,DRYC,SATTH,DRYTH,AA,B1,B2,TMOIST,LOG10E,     &
     &           DHMAX,DHFACT,QHCTOP,OUTTIM,DMAXHE,HTOP,TGRAD,UP,DOWN,  &
     &           STOPHR,GRAV,TSMEAN,TSAMP,DLAI,BARE,HPR
      WRITE(LUB) (Z(I),H(I),NTROOT(I),THETA(I),KL(I),C(I),T(I),         &
     &           CHSOIL(I),I=1,NPT),THETAW,THETAD,THETAN,RDF,FPET,      &
     &           MGR,VAPDIF,VC,TSOIL,PETPC,QWLEAK,QHLEAK
      WRITE(LUB) SWPA,SHPA,SINLAT,COSLAT,TANLAT,TPIY,DAYSEXT,SB,CHW,    &
     &           ALBEDO,ALT,PMB,ZU,ZH,ZM,ZT,D,VK,SEXT,CHA,              &
     &           TCON,CHS,EF,WATDEN,LHV0,CHV,                           &
     &           T0,HSTORE,ISMETH,DHTOL,RHA,IHYS,HYSHPH,SARWA,IPATHA,   &
     &           AIRTOL,HYSTOL,HYSMXH,INDEXB,BBHEAD,HM,THTA,ALFACT
      WRITE(LUB) (PMFN(I),I=1,NYEARS)
      WRITE(LUB) (PRFN(I),I=1,NYEARS)
      WRITE(LUB) PTRANS,PEVAPO
!
!     Read precipitation data; write PET/MET and precip. data to LUB
!
      IF (IFLIST .GT. 0) THEN
        CLOSE(UNIT=LUI)
        IF (IRAIN .EQ. 0) OPEN(UNIT=LUI,FILE=PRFN(1),STATUS='OLD')
      ENDIF
      CALL PRECIPBC(IPRINT,LUSLUW,LUI,LUB)
      CLOSE (UNIT=LUI)
      CLOSE (UNIT=LUB)
      WRITE(LUS,'(80("-"),/,A)') ' Program DATAINH terminated normally.'
      IF (LUSLUW) WRITE(LOUT,'(80("-"),/,A)')                           &
     &' Program DATAINH terminated normally.'
      STOP
999   WRITE(LUS,'(80("-"),/,A)')                                        &
     &' Program DATAINH terminated prematurely.'
      STOP
!-----------------------------------------------------------------------
!     End of program DATAINH
!-----------------------------------------------------------------------
      END
      SUBROUTINE CALPET
!-----------------------------------------------------------------------
!     Called by ETBC
!
!     Calculates PET using Penman formula from FAOPET
!     Note that negative PET values are reported as zero values.
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      REAL*8 NRATIO
      FUNED(TDEW) = EXP(54.878919-(6790.4985/TDEW)-5.02808*LOG(TDEW))
!
!     From Doorenbos and Pruitt (1977), the following equation is 
!     drived from long time mean PMB for a number of stations in Africa
!     (Climate of Africa, Vol 10, World Survey of Climatology)
!
      IF (PMB .LE. 0.0) PMB = 1013.0-0.1152*ALT+5.44*10.E-6*ALT*ALT
      A = 0.17
      IF (ZU .LE. 2.0) A = 0.22
      UHTCF = (2.0/ZU)**A
      SIGMA = 2.0E-9
      FACTRS = 0.017
      FU24 = 38.6
      I = 0
100   I = I+1
      TMAX = RMDATA(2,I)
      TMIN = RMDATA(3,I)
      TDEW = RMDATA(4,I)
      SR_MEAS = RMDATA(5,I)*FACTRS
      WIND  = RMDATA(6,I)*UHTCF*FU24
      CLOUD = RMDATA(7,I)
      TMAX  = (TMAX-32.0)*5.0/9.0
      TMIN  = (TMIN-32.0)*5.0/9.0
      TMEAN = 273.16+(TMAX+TMIN)/2.0
      TDEW  = 273.16+(TDEW-32.0)*5.0/9.0
      EA = FUNED(TMEAN)
      ED = FUNED(TDEW)
      NRATIO = 0.95-0.066*CLOUD-0.0023*CLOUD**2.0
      GG  = 0.0006595*PMB
      D   = (EA/TMEAN)*(6790.4985/TMEAN-5.02808)
      W   = D/(D+GG)
      FT  = SIGMA*TMEAN**4.0
      FED = 0.34-0.044*SQRT(ED)
      FNN = 0.1+0.9*NRATIO
      RN  = (1.0-ALBEDO)*SR_MEAS-FT*FED*FNN
      FU  = 0.27*(1.0+WIND/100.0)
      PET(I) = (W*RN+(1.0-W)*FU*(EA-ED))/10.0
      IF (PET(I) .LT. 0.0) PET(I) = 0.0
      IF (I .LT. NDAYS) GO TO 100
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine CALPET
!-----------------------------------------------------------------------
      END
      SUBROUTINE MYHRLY(FPT)
!-----------------------------------------------------------------------
!     Called by ETPAR
!
!     Calculates hourly distribution of PET based on the sine wave
!     approach of Hillel (SSSAJ 40:807-815, 1976)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION FPT(24)
      TPI = 2.0*PI/24.0
      DO IHOUR=1,24
        IF (IHOUR .LT. 7 .OR. IHOUR .GE. 19) THEN
          FPT(IHOUR) = 0.01
        ELSE
          TEMP = -0.5*COS(TPI*(IHOUR-6))+0.5*COS(TPI*(IHOUR-7))
          FPT(IHOUR) = 0.88*TEMP
        ENDIF
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine MYHRLY
!-----------------------------------------------------------------------
      END
      FUNCTION RLD(DEPTH)
!-----------------------------------------------------------------------
!     Called by PLANTIN and ROOTS
!
!     Calculates root density as a function of depth
!     Note that roots are not allowed in node 1 (the surface node,
!     where DEPTH = 0)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      RLD = AA*EXP(-B1*DEPTH)+B2
      IF (DEPTH .EQ. 0.0) RLD = 0.0
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine RLD
!-----------------------------------------------------------------------
      END
      SUBROUTINE THERMK(IBEGIN,IEND)
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     Calculates soil thermal conductivity and vol. heat capacity
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DO I=IBEGIN,IEND
        MATI  = MAT(I)
        THETAI = THETA(I)
        TCONA = TCON(1,MATI)
        TCONB = TCON(2,MATI)
        TCONC = TCON(3,MATI)
        TCOND = TCON(4,MATI)
        TCONE = TCON(5,MATI)
        KH(I) = TCONA+TCONB*THETAI-TCOND*EXP(-1.0*(TCONC*THETAI)**TCONE)
        CH(I) = CHW*THETAI+CHSOIL(I)
        IF (IVAPOR .EQ. 1) CH(I) = CH(I)+CHV*THETAV(I)
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine THERMK
!-----------------------------------------------------------------------
      END
      SUBROUTINE TIMEX
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     Queries the operating system and returns the time and date
!     SDATE = Date in format 01 Jan 1993 (dd mmm yyyy)
!     STIME = Time in format 01:02:03.42 (hh:mm:ss.ss)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      CHARACTER*3 MON(12)*3,REAL_CLOCK(3)*12
      INTEGER DATE_TIME(8)
      DATA MON /'Jan','Feb','Mar','Apr','May','Jun','Jul',              &
     & 'Aug','Sep','Oct','Nov','Dec'/
      CALL DATE_AND_TIME(REAL_CLOCK(1),REAL_CLOCK(2),                   &
     &  REAL_CLOCK(3),DATE_TIME)
      SDATE(1:2)  = REAL_CLOCK(1)(7:8)
      SDATE(3:3)  = CHAR(32)
      SDATE(4:6)  = MON(DATE_TIME(2))
      SDATE(7:7)  = CHAR(32)
      SDATE(8:11) = REAL_CLOCK(1)(1:4)
      STIME(1:2)  = REAL_CLOCK(2)(1:2)
      STIME(3:3)  = CHAR(58)
      STIME(4:5)  = REAL_CLOCK(2)(3:4)
      STIME(6:6)  = CHAR(58)
      STIME(7:11) = REAL_CLOCK(2)(5:9)
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine TIMEX
!-----------------------------------------------------------------------
      END
      SUBROUTINE POLYKH(IBEGIN,IEND,HEAD,RKL)
!-----------------------------------------------------------------------
!     Called by DATAINH, HDRYCALC, and UNSATH
!
!     Calculates the conductivity (RKL) of each node as a function of 
!     the head value (HEAD) for that node
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION HEAD(*),RKL(*)
      LOGICAL LMDIFF
      MM = 0
      DO I=IBEGIN,IEND
      M = MAT(I)
      IF (M .NE. MM) THEN
        LMDIFF = .TRUE.
        MM = M
      ENDIF
      X = HEAD(I)
      IF (X .LE. AIRINK(M)) THEN 
        RKL(I) = SK(M)
        GO TO 1000
      ENDIF
      GO TO (100,200,300,400,500,600,700,800), KOPT
!-----------------------------------------------------------------------
!     KOPT = 1: Polynomial description of conductivity function
!-----------------------------------------------------------------------
100   IF (LMDIFF) THEN
        IPREV = (M-1)*MAXPOL
        IMAX  = IPREV+NSUBKH(M)
        LMDIFF = .FALSE.
      ENDIF
      J = IPREV+1
110   IF (X .GT. XDIVKH(J) .AND. J .LT. IMAX) THEN
        J = J+1
        GO TO 110
      ENDIF
      IBGN = (M-1)*INC+(J-IPREV-1)*MAXCOE+1
      IND  = IBGN+NDEGKH(J)-1
      XX  = LOG10(X)
      SUM = CREGKH(IND-1)+CREGKH(IND)*XX
      DO K=IND-2,IBGN,-1
        SUM = SUM*XX + CREGKH(K)
      ENDDO
      RKL(I) = 10.0**SUM
      GO TO 1000
!-----------------------------------------------------------------------
!     KOPT = 2: Haverkamp conductivity function 
!-----------------------------------------------------------------------
200   IF (LMDIFF) THEN
        IP = (M-1)*INC
        A  = CREGKH(IP+1)
        B  = CREGKH(IP+2)
        SA = SK(M)*A
        LMDIFF = .FALSE.
      ENDIF
      RKL(I) = SA/(A+X**B)
      GO TO 1000 
!-----------------------------------------------------------------------
!     KOPT = 3: Brooks-Corey w/Burdine or Mualam conductivity functions
!-----------------------------------------------------------------------
300   IF (LMDIFF) THEN
        IP    = (M-1)*INC
        RKMOD = CREGKH(IP+1)
        B     = CREGKH(IP+2)
        EPIT  = CREGKH(IP+3)
        SKM   = SK(M)
        AIRK  = AIRINK(M)
        IF (RKMOD .EQ. 1) THEN
          EXPON = 2.0+(1.0+EPIT)/B
        ELSE
          EXPON = 2.0+(2.0+EPIT)/B
        ENDIF
        LMDIFF = .FALSE.
      ENDIF
      RKL(I) = SKM*(AIRK/X)**EXPON
      GO TO 1000
!-----------------------------------------------------------------------
!     KOPT = 4: van Genuchten w/Burdine or Mualem conductivity functions
!-----------------------------------------------------------------------
400   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        VGA    = CREGKH(IP+1)
        VGN    = CREGKH(IP+2)
        VGM    = CREGKH(IP+3)
        RKMOD  = CREGKH(IP+4)
        EPIT   = CREGKH(IP+5)
        SKM    = SK(M)
        EXPON  = -VGM*EPIT
        LMDIFF = .FALSE.
      ENDIF
      TEMP = 1.0+(VGA*X)**VGN
      IF (RKMOD .EQ. 1) THEN
        RKL(I) = SKM*TEMP**EXPON*(1.0-(1.0-1.0/TEMP)**VGM)
      ELSE
        RKL(I) = SKM*TEMP**EXPON*(1.0-(1.0-1.0/TEMP)**VGM)**2
      ENDIF
      GO TO 1000
!-----------------------------------------------------------------------
!     KOPT = 5: Fayer-Simmons Brooks-Corey w/Mualem function
!-----------------------------------------------------------------------
500   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        B      = CREGTH(IP+1)
        AIRTM  = AIRINT(M)
        THETM  = THET(M)
        THTAM  = THTA(M)
        EXPO   = -1.0/B
        HMM    = HM(M)
        RLOGHM = LOG(HMM)
        THTAHM = THTAM/RLOGHM
        EPIT   = CREGKH(IP+1)
        BC1    = CREGKH(IP+2)
        BC2    = CREGKH(IP+3)
        BC3    = CREGKH(IP+4)
        PM1    = CREGKH(IP+5)
        GAMMAM = CREGKH(IP+6)
        GAMMAX = CREGKH(IP+7)
        SKM    = SK(M)
        AIRK   = AIRINK(M)
        LMDIFF = .FALSE.
      ENDIF
      TEMP   = X/AIRTM
      TEMPE  = TEMP**EXPO
      CHITA  = (1-LOG(X)/RLOGHM)*THTAM
      RTHETA = CHITA+(THETM-CHITA)*TEMPE
      S      = RTHETA/THETM
      GAMMAS = (BC3+BC2*LOG(X))*X**PM1/PM1 + BC1/X - GAMMAM
      RKL(I) = SKM*S**EPIT*(GAMMAS/GAMMAX)**2
      GO TO 1000
!-----------------------------------------------------------------------
!     KOPT = 6: Fayer-Simmons van Genuchten w/Mualem function
!-----------------------------------------------------------------------
600   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        VGA    = CREGTH(IP+1)
        VGN    = CREGTH(IP+2)
        VGM    = CREGTH(IP+3)
        EPIT   = CREGKH(IP+1)
        WO     = CREGKH(IP+2)
        GAMMAM = CREGKH(IP+3)
        GAMMA  = CREGKH(IP+4)
        RLAMB  = CREGKH(IP+5)
        SA     = CREGKH(IP+6)
        WHM    = CREGKH(IP+7)
        CK1    = CREGKH(IP+8)
        CK2    = CREGKH(IP+9)
        CK3    = CREGKH(IP+10)
        HXO    = CREGKH(IP+11)
        GSXO   = CREGKH(IP+12)
        RI3BWO = CREGKH(IP+13)
        GAMMAX = CREGKH(IP+14)
        SKM    = SK(M)
        THETM  = THET(M)
        THTAM  = THTA(M)
        HMM    = HM(M)
        SA     = THTAM/THETM
        P      = 1+VGM
        LMDIFF = .FALSE.
      ENDIF
      W   = 1.0/(1.0+(VGA*X)**VGN)
      XSA = SA-LOG(X)*SA*RLAMB
      S   = XSA+(1-XSA)*W**VGM
      CALL INTGRL(X,W,GAMMAS)
      RKL(I) = SKM*S**EPIT*(GAMMAS/GAMMAX)**2
      GO TO 1000
!-----------------------------------------------------------------------
!     KOPT = 7: Rossi-Nimmo "sum" w/ Mualem conductivity
!-----------------------------------------------------------------------
700   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        EPIT   = CREGKH(IP+1)
        G2CON  = CREGKH(IP+2)
        GAMMA2 = CREGKH(IP+3)
        GAMMAX = CREGKH(IP+4)
        RNA    = CREGTH(IP+1)
        RNC    = CREGTH(IP+2)
        PSID   = CREGTH(IP+3)
        PSIO   = CREGTH(IP+4)
        RLAM   = CREGTH(IP+5)
        PSII   = CREGTH(IP+6)
        SKM    = SK(M)
        LMDIFF = .FALSE.
      ENDIF
      IF (X .GE. PSII) THEN
        TOTAL = ((RNA+(RLAM/(1+RLAM))*(PSIO/X)**RLAM)/X)-G2CON
      ELSE
        TOTAL = GAMMA2+2*RNC*(PSII-X)/PSIO**2
      ENDIF
      RKL(I) = SKM*(THETA(I)/THET(M))**EPIT*(TOTAL/GAMMAX)**2
      GO TO 1000
!-----------------------------------------------------------------------
!     KOPT = 8: Rossi-Nimmo "junction" w/Burdine or Mualem conductivity
!-----------------------------------------------------------------------
800   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        EPIT   = CREGKH(IP+1)
        BETA   = CREGKH(IP+2)
        GAMMA3 = CREGKH(IP+3)
        GAMMA2 = CREGKH(IP+4)
        GAMMAX = CREGKH(IP+5)
        RNA    = CREGTH(IP+1)
        RNC    = CREGTH(IP+2)
        PSID   = CREGTH(IP+3)
        PSIO   = CREGTH(IP+4)
        RLAM   = CREGTH(IP+5)
        PSII   = CREGTH(IP+6)
        PSIJ   = CREGTH(IP+7)
        SKM    = SK(M)
        LMDIFF = .FALSE.
      ENDIF
      IF (X .GE. PSIJ) THEN
        TOTAL = RNA*((1.0/X)-1.0/PSID)
      ELSE IF (X .GE. PSII) THEN
        TOTAL = GAMMA3+BETA*(X**(-RLAM-1)-PSIJ**(-RLAM-1))
      ELSE
        TOTAL = GAMMA3+GAMMA2+2*RNC*(PSII-X)/PSIO**2
      ENDIF
      RKL(I) = SKM*(THETA(I)/THET(M))**EPIT*(TOTAL/GAMMAX)**2
!
!     End of conductivity options
!
1000    IF (RKL(I) .LE. 0.0) THEN
          WRITE(*,6000) IN,I,IDAY,TSTEPB,HH(I),H(I),RKL(I)
          STOP
        ENDIF
      ENDDO
      RETURN
6000  FORMAT(/,
     &' Subroutine POLYKH detected an invalid conductivity value',/,    &
     &I13,' = Iteration number (IN)',/,                                 &
     &I13,' = Node number (N)',/,                                       &
     &I13,' = Simulation day (IDAY)',/,                                 &
     &1PE13.5,' = Completed time of day (TSTEPB)',/,                    &
     &E13.5,' = Previous head HH(I)',/,                                 &
     &E13.5,' = Current head H(I)',/,                                   &
     &E13.5,' = Calculated conductivity RKL(I)',/)
!-----------------------------------------------------------------------
!     End of subroutine POLYKH
!-----------------------------------------------------------------------
      END
      SUBROUTINE RETENT(IBEGIN,IEND,HEAD,IFLAGT,RTHETA,IFLAGC,RC)
!-----------------------------------------------------------------------
!     Called by DATAINH, PLANTIN, UNSATH, and HDRYCALC
!
!     Calculates the water content (THETA) and moisture capacity (C) 
!     of each node as a function of the head value (H) for that node
!       IFLAGC = flag signifying that capacity will be calculated
!       IFLAGT = flag signifying that water content will be calculated
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION HEAD(*),RTHETA(*),RC(*)
      LOGICAL LMDIFF
      INCLUDE 'init.inc'
      MM = 0
      DO 3000 I=IBEGIN,IEND
      M = MAT(I)
      IF (M .NE. MM) THEN
        LMDIFF = .TRUE.
        MM = M
      ENDIF
      X = HEAD(I)
      IF (X .LE. AIRINT(M)) THEN
        IF (IFLAGC .EQ. 1) RC(I) = -WCOMPR*THET(M)
        IF (IFLAGT .EQ. 1) RTHETA(I) = THET(M)              
        GO TO 3000
      ENDIF
      GO TO (100,200,300,400,500,600,700,800), KOPT
!-----------------------------------------------------------------------
!     KOPT = 1: Polynomial description of retention/capacity function
!-----------------------------------------------------------------------
100   IF (LMDIFF) THEN
        IPREV = (M-1)*MAXPOL
        IMAX  = IPREV+NSUBTH(M)
        LMDIFF = .FALSE.
      ENDIF
      J = IPREV+1
110   IF (X .GT. XDIVTH(J) .AND. J .LT. IMAX) THEN
        J = J+1
        GO TO 110
      ENDIF
      IBGN = (M-1)*INC+(J-IPREV-1)*MAXCOE+1
      IND  = IBGN+NDEGTH(J)-1
      XX   = LOG10(X)
      IF (IFLAGT .EQ. 1) THEN
        SUM = CREGTH(IND-1)+CREGTH(IND)*XX
        DO K=IND-2,IBGN,-1
          SUM = SUM*XX+CREGTH(K)
        ENDDO
        RTHETA(I) = SUM
      ENDIF
      IF (IFLAGC .EQ. 1) THEN
        SUM = (IND-1-IBGN)*CREGTH(IND-1)+(IND-IBGN)*CREGTH(IND)*XX
        DO K=IND-2,IBGN+1,-1
          SUM = SUM*XX+(K-IBGN)*CREGTH(K)
        ENDDO
        RC(I)  = LOG10E*SUM/X
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 2: Haverkamp retention/capacity function
!-----------------------------------------------------------------------
200   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        ALPHA  = CREGTH(IP+1)
        BETA   = CREGTH(IP+2)
        RETOPT = CREGTH(IP+3)
        TR     = THTR(M)
        THETR  = THET(M)-TR
        LMDIFF = .FALSE.
      ENDIF
      IF (RETOPT .EQ. 1) THEN
        XX  = X
      ELSE
        XX  = LOG(X)
      ENDIF
      XBETA = XX**BETA
      TEMP  = ALPHA*THETR/(ALPHA+XBETA)
      IF (IFLAGT .EQ. 1) RTHETA(I) = TEMP+TR
      IF (IFLAGC .EQ. 1) THEN
        IF (RETOPT .EQ. 1) THEN
          RC(I) = -(XX**(BETA-1))*ALPHA*BETA*THETR/(ALPHA+XBETA)**2
        ELSE IF (RETOPT .EQ. 2) THEN
          RC(I)= -(XX**(BETA-1))*ALPHA*BETA*THETR/(X*(ALPHA+XBETA)**2)
        ENDIF
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 3: Brooks-Corey retention/capacity function
!-----------------------------------------------------------------------
300   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        B      = CREGTH(IP+1)
        TR     = THTR(M)
        THETR  = THET(M)-TR
        AIRT   = AIRINT(M)
        EXP    = -1.0/B
        LMDIFF = .FALSE.
      ENDIF
      TEMP = THETR*(X/AIRT)**EXP
      IF (IFLAGT .EQ. 1) RTHETA(I) = TEMP+TR
      IF (IFLAGC .EQ. 1) THEN
        RC(I) = -TEMP/(B*X)
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 4: van Genuchten retention/capacity function 
!-----------------------------------------------------------------------
400   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        VGA    = CREGTH(IP+1)
        VGN    = CREGTH(IP+2)
        VGM    = CREGTH(IP+3)
        AMN    = VGA*VGM*VGN
        THETM  = THET(M)
        THTRM  = THTR(M)
        THETR  = THETM-THTRM
        LMDIFF = .FALSE.
      ENDIF
      TEMP  = 1.0+(X*VGA)**VGN
      TEMPM = TEMP**(-VGM)
      IF (IFLAGT .EQ. 1) RTHETA(I) = THTRM+THETR*TEMPM
      IF (IFLAGC .EQ. 1) THEN
        AHNM1 = (VGA*X)**(VGN-1)
        RC(I) = -THETR*AMN*AHNM1*TEMPM/TEMP
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 5: Fayer-Simmons Brooks-Corey retention/capacity function
!-----------------------------------------------------------------------
500   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        B      = CREGTH(IP+1)
        AIRTM  = AIRINT(M)
        THETM  = THET(M)
        THTAM  = THTA(M)
        EXP    = -1.0/B
        HMM    = HM(M)
        RLOGHM = LOG(HMM)
        THTAHM = THTAM/RLOGHM
        LMDIFF = .FALSE.
      ENDIF
      TEMP  = X/AIRTM
      TEMPE = TEMP**EXP
      CHITA = (1-LOG(X)/RLOGHM)*THTAM
      IF (IFLAGT .EQ. 1) RTHETA(I) = CHITA+(THETM-CHITA)*TEMPE
      IF (IFLAGC .EQ. 1) THEN
        OOBHE = 1/(B*AIRTM)
        RC(I) = (TEMPE-1)*THTAHM/X - (THETM-CHITA)*OOBHE*TEMPE/TEMP
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 6: Fayer-Simmons van Genuchten retention/capacity function 
!-----------------------------------------------------------------------
600   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        VGA    = CREGTH(IP+1)
        VGN    = CREGTH(IP+2)
        VGM    = CREGTH(IP+3)
        AMN    = CREGTH(IP+6)
        RLOGHM = CREGTH(IP+7)
        THTAHM = CREGTH(IP+8)
        THETM  = THET(M)
        THTAM  = THTA(M)
        LMDIFF = .FALSE.
      ENDIF
      TEMP  = 1.0+(X*VGA)**VGN
      TEMPM = TEMP**(-VGM)
      CHITA = (1-LOG(X)/RLOGHM)*THTAM
      IF (IFLAGT .EQ. 1) RTHETA(I) = CHITA+(THETM-CHITA)*TEMPM
      IF (IFLAGC .EQ. 1) THEN
        AHNM1 = (VGA*X)**(VGN-1)
        RC(I) = (TEMPM-1)*THTAHM/X-(THETM-CHITA)*AMN*AHNM1*TEMPM/TEMP
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 7: Rossi-Nimmo "sum" retention/capacity function 
!-----------------------------------------------------------------------
700   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        RNA    = CREGTH(IP+1)
        RNC    = CREGTH(IP+2)
        PSID   = CREGTH(IP+3)
        PSIO   = CREGTH(IP+4)
        RLAM   = CREGTH(IP+5)
        PSII   = CREGTH(IP+6)
        THETM  = THET(M)
        LMDIFF = .FALSE.
      ENDIF
      IF (IFLAGT .EQ. 1) THEN
        IF (X .LE. PSII) THEN
          RTHETA(I) = THETM*(1.0-RNC*(X/PSIO)**2)
        ELSE
          RTHETA(I) = THETM*((PSIO/X)**RLAM-(PSIO/PSID)**RLAM           &
     &                +RNA*LOG(PSID/X))
        ENDIF
      ENDIF
      IF (IFLAGC .EQ. 1) THEN
        IF (X .LE. PSII) THEN
          RC(I) = -THETM*2.0*RNC*X/PSIO**2
        ELSE
          RC(I) = -(THETM/X)*(RNA+RLAM*(PSIO/X)**RLAM)
        ENDIF
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 8: Rossi-Nimmo "junction" retention/capacity function 
!-----------------------------------------------------------------------
800   IF (LMDIFF) THEN
        IP     = (M-1)*INC
        RNA    = CREGTH(IP+1)
        RNC    = CREGTH(IP+2)
        PSID   = CREGTH(IP+3)
        PSIO   = CREGTH(IP+4)
        RLAM   = CREGTH(IP+5)
        PSII   = CREGTH(IP+6)
        PSIJ   = CREGTH(IP+7)
        THETM  = THET(M)
        LMDIFF = .FALSE.
      ENDIF
      IF (IFLAGT .EQ. 1) THEN
        IF (X .LE. PSII) THEN
          RTHETA(I) = THETM*(1.0-RNC*(X/PSIO)**2)
        ELSEIF (X .GE. PSIJ) THEN
          RTHETA(I) = THETM*RNA*LOG(PSID/X)
        ELSE
          RTHETA(I) = THETM*(PSIO/X)**RLAM
        ENDIF
      ENDIF
      IF (IFLAGC .EQ. 1) THEN
        IF (X .LE. PSII) THEN
          RC(I) = -THETM*2.0*RNC*X/PSIO**2
        ELSEIF (X .GE. PSIJ) THEN
          RC(I) = -THETM*RNA/X
        ELSE
          RC(I) = -THETM*RLAM*((PSIO/X)**RLAM)/X
        ENDIF
        GO TO 2000
      ELSE
        GO TO 3000
      ENDIF
!
!     Capacity Error Message
!
2000  IF (RC(I) .GT. 0.0) THEN
        IERROR = IERROR + 1
        WRITE(*,6000) M,I,-RC(I),HEAD(I)
      ENDIF
3000  CONTINUE
      IF (IFLAGC .EQ. 1) THEN
      IF (IERROR .GT. 0) THEN
        WRITE(*,6010) IERROR
        STOP
      ENDIF
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
6000  FORMAT(                                                           &
     &' ERROR: Unsaturated capacity term negative, MAT No.',I3,/,       &
     &'        Node = ',I4,', RC = ',1PE15.7,', HEAD = ',E15.7)
6010  FORMAT(/,                                                         &
     &' Program execution stopped in subroutine RETENT ',/,             &
     &' because of negative capacity term(s) for',I4,' nodes')
!-----------------------------------------------------------------------
!     End of Subroutine RETENT
!-----------------------------------------------------------------------
      END
!-----------------------------------------------------------------------
!     File KCALCS
!
!     This file contains the subroutines INTGRL, FGX, FHX, and FGSX
!
!     Subroutines are called as follows:
!        INTGRL called by POLYKH and SHPPAR
!        FGX called by INTGRL
!        FHX called by FGX and SHPPAR
!        FGSX called by FGX and SHPPAR
!
!     Calculate hydraulic conductivity for the Fayer-Simmons retention
!     function. Some of the preliminary calculations are performed
!     in DATAINH. See Fayer and Simmons (1995) 31:1233-1238.
!-----------------------------------------------------------------------
      SUBROUTINE INTGRL(HEAD,W,GAMMAS)
!-----------------------------------------------------------------------
!     Calculates GAMMAS in integration of new retention function
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      CALL FGX(W,G1W)
      CALL FGX(WO,G1WO)
      RI1 = (1.-WHM)**VGM - (1.-W)**VGM
      RI2 = (1./HEAD-1./HMM)/VGA                                        &
     &      + (1.-WHM)**(VGM-1.) - (1.-W)**(VGM-1.)
      IF (W .LE. WO) THEN
        RI3B = ((W-WHM)+WHM*LOG(WHM)-W*LOG(W))*VGM
        RI3  = RI3B
      ELSE
        RI3A = (1.-WO)**VGM*(LOG((1.-WO)/WO)-1./VGM)                    &
     &         -(1.-W)**VGM*(LOG((1.-W)/W)-1./VGM) + G1W - G1WO
        RI3  = RI3A+RI3BWO
      ENDIF
      GAMMAS = RI1*VGA*(1.-GAMMA*SA-SA)                                 &
     &         + (VGA*RLAMB*SA)*(RI2+RI3/VGN) + GAMMAM
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine INTGRL
!-----------------------------------------------------------------------
      END
      SUBROUTINE FGX(X,GX)
!-----------------------------------------------------------------------
!     Calculates G function in integration of new retention function
!     X is either w or wo
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (1.-X .GT. XO) THEN
        CALL FHX(X,HX)
        GX = GSXO + HXO - HX
      ELSE
        CALL FGSX(1.-X,GSX)
        GX = GSX
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine FGX
!-----------------------------------------------------------------------
      END
      SUBROUTINE FHX(X,HX)
!-----------------------------------------------------------------------
!     Calculates H function in integration of new retention function
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      X2 = X*X
      HX = LOG(X)+CK1*X+CK2*X2/2+CK3*X2*X/3
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine FHX
!-----------------------------------------------------------------------
      END
      SUBROUTINE FGSX(X,GSX)
!-----------------------------------------------------------------------
!     Calculates GSX function in integration of new retention function
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      SUM = X/P
      PSUM = SUM
      DO K=1,100
        SUM = SUM+X**K/(P+K)
        IF (SUM/PSUM .LT. 1.0001) GO TO 100
        PSUM = SUM
      ENDDO
100   GSX = SUM*X**P
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine FGSX
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     End of file KCALCS (subroutines INTGRL, FGX, FHX, and FGSX)
!-----------------------------------------------------------------------
      END

      SUBROUTINE WELCOME(name,ver)
!-----------------------------------------------------------------------
!     Called by DATAINH, UNSATH, and DATAOUT
!
!     Writes program header with version number and contact person
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      CHARACTER name*7
      WRITE(LUS,10) name,ver
10    FORMAT(/,                                                         &
     &'         Program ',A7,/,                                         &
     &'         Version ',F4.2,//,                                      &
     &'         Contact:',/,                                            &
     &'           MJ Fayer',/,                                          &
     &'           Box 999, MSIN K9-33',/,                               &
     &'           Richland, WA  99352',/,                               &
     &'           phone 509-372-6045',/,                                &
     &'           FAX   509-372-6089',/,                                &
     &'           email mike.fayer@pnl.gov',/)
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine WELCOME
!-----------------------------------------------------------------------
      END


      SUBROUTINE RELHUM(HEAD,TEMP,RHOUT)
!-----------------------------------------------------------------------
!     Called by KVCALC, UNSATH, and VOLVAP
!
!     Calculates the relative humidity using suction head,
!     temperature, and osmotic potential (in suction head units)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (HEAD .LE. 0.0D0) THEN
        RHOUT = 1.0D0
      ELSE
        RHOUT = EXP(-(HEAD+OSMPOT)*MGR/TEMP)
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine RELHUM
!-----------------------------------------------------------------------
      END
      SUBROUTINE VOLVAP(IBEGIN,IEND)
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     Calculates volumetric vapor content (THETAV)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      FUNVD(TD) = EXP(46.440973-(6790.4985/TD)-6.02808*LOG(TD))
      DO I=IBEGIN,IEND
        SVD = FUNVD(T(I))
        CALL RELHUM(H(I),T(I),RHUMID)
        RHOV = SVD*RHUMID
        THETAV(I) = RHOV*(THET(MAT(I))-THETA(I))/WATDEN
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine VOLVAP
!-----------------------------------------------------------------------
      END
      SUBROUTINE HYSINI
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     This subroutine initializes variables for the hysteresis code.
!     It is only called once, at the beginning of a simulation.
!     This code is based on the air-water portion of the Lenhard
!     and Parker air-oil-water hysteresis code.  Water retention is
!     described with the van Genuchten function (assuming m = 1-1/n)
!     and the Mualem conductivity model (assuming a pore interaction 
!     exponent of 0.5)
!
!     ASW    = Effective apparent water saturation
!     ALPHAD = van Genuchten parameter alpha for primary drainage path
!     ALPHAI = van Genuchten parameter alpha for primary wetting path
!     ASAT   = Entrapped air saturation
!     CWW    = Capacity term for water saturation
!     ESW    = Effective water saturation
!     HYSMXH = Head above which hysteresis and air entrapment are neglected
!              Use of HYSMXH is valid as long as HDRY is greater and the
!              options ISHOPT and IHEAT are not set equal to 1
!     HAW    = Air-water capillary head
!     I      = Node
!     IHYS   = Index specifying initial path
!                0 = no hysteresis
!                1 = primary drainage path
!                2 = wetting from primary drainage path
!                3 = drying from wetting path (see IHYS=2)
!                4 = restart from previous simulation using *.HRI file
!     IPSW   = Index for saturation path (i.e., primary drainage = 1,
!              primary wetting = 2, secondary drainage = 3, etc.);
!     IPATH  = Maximum number of wetting and drainage paths (default is
!              seven). IPATH = 1 is the primary drainage path
!     JJH    = Index for drying (0) or wetting (1) path
!     LUH    = Restart data file
!     PERMW  = Relative permeability of water
!     PHI    = Soil porosity (equivalent to saturated water content)
!     RASW   = Reversal apparent water saturation
!     RHSW   = Reversal capillary head associated with RASW
!     RSW    = Effective saturation at primary drainage branch reversal
!     RSWAW  = Effective saturation at primary drainage branch reversal.
!              Determines air entrapment by air-water interface
!     SARW   = Maximum effective entrapped air on a given scanning path
!     SARWA  = Maximum effective entrapped air on the primary imbibition
!              path; specified for each material.
!     SATW   = Effective saturation of air trapped by water
!     SM     = van Genuchten parameter ThetaR expressed as saturation
!     SRAW   = Fraction of apparent water sat. above reversal water sat.
!              that is entrapped water
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      LOGICAL EXT
      OPTOL3 = 1+HYSTOL
      OMTOL3 = 1-HYSTOL
!-----------------------------------------------------------------------
!     Calculate RAW value for each node based on SARWA value for
!     material type of that node.  Also calculate ASW value (ASWMXH)
!     associated with HYSMXH
!-----------------------------------------------------------------------
      DO I=1,NPT
        M = MAT(I)
        RAWA(M) = 0.0
        IF (SARWA(M) .GT. 0) RAWA(M) = -1.0+1.0/SARWA(M)
        IP     = (M-1)*INC
        ALPHAD = CREGTH(IP+1)
        XN     = CREGTH(IP+2)
        XM     = CREGTH(IP+3)
        ASWMXH(I) = (1.0+(ALPHAD*HYSMXH)**XN)**(-XM)
      ENDDO
      IF (IHYS .LT. 4) THEN
!-----------------------------------------------------------------------
!     Initialize reversal points for paths 3 to IPATH
!-----------------------------------------------------------------------
        DO I=1,NPT
          RASW(I,1) = 1.0
          RHSW(I,1) = 0.0
          IPATH = IPATHA(MAT(I))
          DO J=3,IPATH,2
            RASW(I,J) = ASWMXH(I)
            RHSW(I,J) = HYSMXH
          ENDDO
          DO J=2,IPATH,2
            RASW(I,J) = 1.0
            RHSW(I,J) = 0.0
          ENDDO
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     Set initial conditions and reversal points for paths 1
!     to 3 based on starting saturation path.
!-----------------------------------------------------------------------
      IF (IHYS .EQ. 1) THEN
!-----------------------------------------------------------------------
!       Start on the primary drainage path
!-----------------------------------------------------------------------
        DO I = 1,NPT
          JJH(I)   = 0
          IPSW(I)  = 1
          RSWAW(I) = 1.0
        ENDDO
      ELSE IF (IHYS .EQ. 2) THEN
!-----------------------------------------------------------------------
!     Start on imbibition path defined by HYSHPH
!     When HYSHPH equals HYSMXH, the path is the primary imbibition path
!-----------------------------------------------------------------------
        DO I = 1,NPT
          M = MAT(I)
          IP     = (M-1)*INC
          ALPHAD = CREGTH(IP+1)
          XN     = CREGTH(IP+2)
          XM     = CREGTH(IP+3)
          ASWHPH = (1.0+(ALPHAD*HYSHPH(M))**XN)**(-XM)
          RHSW(I,2) = HYSHPH(M)
          RASW(I,2) = ASWHPH
          RSWAW(I) = ASWHPH
          SARW(I)  = (1.0-ASWHPH)/(1.0+RAWA(M)*(1.0-ASWHPH))
          JJH(I)   = 1
          IPSW(I)  = 2
        ENDDO
      ELSE IF (IHYS .EQ. 3) THEN
!-----------------------------------------------------------------------
!     Start on drainage path defined by imbibition path from HYSHPH
!-----------------------------------------------------------------------
        DO I = 1,NPT
          M = MAT(I)
          IP     = (M-1)*INC
          ALPHAD = CREGTH(IP+1)
          XN     = CREGTH(IP+2)
          XM     = CREGTH(IP+3)
          ASWHPH = (1.0+(ALPHAD*HYSHPH(M))**XN)**(-XM)
          RHSW(I,2) = HYSHPH(M)
          RASW(I,2) = ASWHPH
          RHSW(I,3) = 0.0
          RASW(I,3) = 1.0
          RSWAW(I) = ASWHPH
          SARW(I)  = (1.0-ASWHPH)/(1.0+RAWA(M)*(1.0-ASWHPH))
          JJH(I)   = 0
          IPSW(I)  = 3
        ENDDO
      ELSE IF (IHYS .EQ. 4) THEN
!-----------------------------------------------------------------------
!     Restart - read initial and reversal values from restart file
!-----------------------------------------------------------------------
        INQUIRE(FILE=HYFILE,EXIST=EXT)
        IF (EXT) THEN
          OPEN(UNIT=LUH,FILE=HYFILE,STATUS='OLD')
        ELSE
          WRITE(LUS,*) ' The hysteresis restart file was not present'
          GO TO 999
        ENDIF
        READ(LUH,2005)( IPSW(I),I=1,NPT)
        READ(LUH,2005)  (JJH(I),I=1,NPT)
        READ(LUH,2010) (SARW(I),I=1,NPT)
        READ(LUH,2010)(RSWAW(I),I=1,NPT)
        DO I=1,NPT
          IPATH = IPATHA(MAT(I))
          READ(LUH,2010)(RASW(I,J),J=1,IPATH)
          READ(LUH,2010)(RHSW(I,J),J=1,IPATH)
          IF (RHSW(I,2) .GT. HYSMXH) THEN
            WRITE(*,2020) I,RHSW(I,2),HYSMXH
            GO TO 999
          ENDIF
        ENDDO
      ENDIF
      RETURN
999   WRITE(LUS,*) ' Subroutine HYSINI was terminated prematurely'
      STOP
!-----------------------------------------------------------------------
!     FORMAT Statements
!-----------------------------------------------------------------------
2005  FORMAT(1X,70I1)
2010  FORMAT(1P7E15.7)
2020  FORMAT(/,                                                         &
     &' INPUT ERROR when reading the restart file.  Details are:',/,    &
     &' Node = ',I5,', RHSW(I,2) = ',E15.7,', and HYSMXH = ',E15.7)
!-----------------------------------------------------------------------
!     End of subroutine HYSINI
!-----------------------------------------------------------------------
      END
      SUBROUTINE HYSSHP(IBEGIN,IEND,IOFLAG,RKL,RC,RTHETA)
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     Calls PATH and UPDATE
!
!     This subroutine passes the UNSAT-H parameters into the hysteresis
!     routines in the correct units and converts the output to UNSAT-H
!     units.  
!     IBEGIN = Beginning node of calculation loop
!     IEND   = Ending node of calculation loop
!     IOFLAG = Output flag.  If value is 1, call subroutine UPDATE,
!              which updates hysteresis parameters following successful
!              completion of a time step.
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION RKL(*),RC(*),RTHETA(*)
      DO I=IBEGIN,IEND
        M      = MAT(I)
        IP     = (M-1)*INC
        ALPHAD = CREGTH(IP+1)
        XN     = CREGTH(IP+2)
        XM     = CREGTH(IP+3)
        PHI    = THET(M)
        SM     = THTR(M)/PHI
        XXM    = 1.0/XM
        IPATH  = IPATHA(M)
        IF (IPATH .EQ. 1) THEN
          ALPHAI = ALPHAD
        ELSE
          ALPHAI = ALFACT(M)*ALPHAD
        ENDIF
        HAW = H(I)
        IF (HAW .LT. 0.0) HAW = 0.0
        PHAW = HH(I)
        IF (PHAW .LT. 0.0) PHAW = 0.0
        CALL PATH(HAW,PHAW,I)
        RKL(I) = PERMW*SK(M)
        RC(I)  = -CWW
      IF (HAW .LE. AIRINT(M)) RC(I) = -WCOMPR*THET(M)
        RTHETA(I) = SW*THET(M)
        ENTAIR(I) = ASAT*THET(M)
        JH(I) = JJJH
        IF (IOFLAG .EQ. 1) CALL UPDATE(HAW,PHAW,I)
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine HYSSHP
!-----------------------------------------------------------------------
      END
      SUBROUTINE PATH(HAW,PHAW,I)
!-----------------------------------------------------------------------
!     Called by HYSSHP
!
!     Calculates the air-water capillary head, the amount of entrapped
!     air, and the water saturation using the input pressure heads of
!     water and air
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
!
!     Transfer global variables to local variables.
!
      JJJH  = JJH(I)
      IPSW1 = IPSW(I)
      RAW   = RAWA(MAT(I))
      RSW   = RSWAW(I)
      SRAW  = SARW(I)
      LSW   = 0
!
!     Check the value of IPATH to enter proper calculation loop
!     IPATH = 1 is the primary drainage path (no hysteresis)
!
      IF (IPATH .EQ. 1) THEN
        CALL DRAIN2P(HAW,PHAW,I)
        RETURN
      ENDIF
!
!     Hysteresis Active
!
!     Determine path, then call DRAIN2P, DRY2P, or WET2P
!
!     Check to see if on primary drainage path
!
      IF (HAW .GE. RHSW(I,2)) THEN
        CALL DRAIN2P(HAW,PHAW,I)
        JJJH  = 0
        IPSW1 = 1
        RETURN
      ENDIF
!
!     Check to see if on wetting path
!       a) wetting from a drying path
!       b) wetting on a wetting path
!       c) wetting along the maximum wetting path
!       d) saturated along a wetting path
!
60    IF ((((HAW .LT. PHAW*OMTOL3 .AND. JJJH .EQ. 0)                    &
     &  .OR.(HAW .LT. PHAW*OPTOL3 .AND. JJJH .EQ. 1))                   &
     &.AND. IPSW1 .LT. IPATH) .OR. (JJJH .EQ. 1                         & 
     &.AND. (HAW .EQ. 0.0 .OR. IPSW1 .EQ. IPATH)))THEN
!
!       if on maximum path and node dries above RHSW(I,IPATH),
!       reset ipsw1 and jjjh and send to dry path calculations
!
        IF (IPSW1 .EQ. IPATH .AND. HAW .GE. RHSW(I,IPATH)) THEN
          IPSW1 = IPATH-1
          JJJH  = 0
          GO TO 70
        ENDIF
!
!       check for path reversal from dry to wet; update IPSW
!
        IF (JJJH .EQ. 0 .AND. IPSW1 .LT. IPATH) THEN
          IPSW1 = IPSW1+1
          JJJH  = 1
        ENDIF
!
!       define applicable path and call WET2P; set LSW and JJJH
!
        DO L = 0,IPSW1-1,2
          IF (HAW.LE.RHSW(I,IPSW1-L).AND.HAW.GE.RHSW(I,IPSW1-1-L)) THEN
            IPSW1 = IPSW1-L
            CALL WET2P(HAW,PHAW,I)
            LSW  = L
            JJJH = 1
            RETURN
          ENDIF
        ENDDO
      ENDIF
!
!     Check to see if on drying path
!       a) drying from a wetting path
!       b) drying on a drying path
!       c) drying along the maximum drying path
!
70    IF ((((HAW .GE. PHAW*OPTOL3 .AND. JJJH .EQ. 1)                    &
     & .OR. (HAW .GE. PHAW*OMTOL3 .AND. JJJH .EQ. 0))                   &
     & .AND. IPSW1.LT.IPATH)                                            &
     & .OR. (JJJH .EQ. 0 .AND. IPSW1 .EQ. IPATH)) THEN
!
!       check whether on primary drainage path
!
        IF (IPSW1 .EQ. 1) THEN
          CALL DRAIN2P(HAW,PHAW,I)
          JJJH  = 0
          IPSW1 = 1
          RETURN
        ENDIF
!
!       if on maximum path and node wets below RHSW(I,IPATH),
!       reset ipsw1 and jjjh and send to wet path calculations
!
        IF (IPSW1 .EQ. IPATH .AND. HAW .LE. RHSW(I,IPATH)) THEN
          IPSW1 = IPATH-1
          JJJH  = 1
          GO TO 60
        ENDIF
!
!       check for path reversal from wet to dry; update IPSW
!
        IF (JJJH .EQ. 1 .AND. IPSW1 .LT. IPATH) THEN
          IPSW1 = IPSW1+1
          JJJH  = 0
        ENDIF
!
!       define applicable loop and call DRY2P; set LSW and JJJH
!
        DO L = 0,IPSW1-2,2
          IF(HAW.GE.RHSW(I,IPSW1-L).AND.HAW.LE.RHSW(I,IPSW1-1-L))THEN
            IPSW1 = IPSW1-L
            CALL DRY2P(HAW,PHAW,I)
            LSW  = L
            JJJH = 0
            RETURN
          ENDIF
        ENDDO
      ENDIF
!
!       The IF statements were somehow bypassed
!
      WRITE(LUS,2000) IDAY,HOUR+TSUB-DELT,I,IPSW(I),IPSW1,JJH(I),JJJH,  &
     &  HAW,PHAW,HAW-PHAW,                                              &
     &  RASW(I,1),RHSW(I,1),RASW(I,2),RHSW(I,2),RASW(I,3),RHSW(I,3),    &
     &  RASW(I,4),RHSW(I,4),RASW(I,5),RHSW(I,5),RASW(I,6),RHSW(I,6),    &
     &  RASW(I,7),RHSW(I,7)
      STOP
!-----------------------------------------------------------------------
!     FORMAT Statements
!-----------------------------------------------------------------------
2000  FORMAT(' The IF statements in PATH were somehow bypassed',/,      &
     &' IDAY = ',I3,' TDONE = ',1PE12.5,/,                              &
     &' I    = ',I3,'  IPSW = ',I2,'  IPSW1 = ',I2,/,                   &
     &' JJH  = ',I3,'  JJJH = ',I2,/,                                   &
     &' HAW   = ',1PE13.6,'  PHAW  = ',E13.6,'  HAW-PHAW = ',E13.6,/,   &
     &' RASW1 = ',E13.6,10X,'RHSW1 = ',E13.6,/,                         &
     &' RASW2 = ',E13.6,10X,'RHSW2 = ',E13.6,/,                         &
     &' RASW3 = ',E13.6,10X,'RHSW3 = ',E13.6,/,                         &
     &' RASW4 = ',E13.6,10X,'RHSW4 = ',E13.6,/,                         &
     &' RASW5 = ',E13.6,10X,'RHSW5 = ',E13.6,/,                         &
     &' RASW6 = ',E13.6,10X,'RHSW6 = ',E13.6,/,                         &
     &' RASW7 = ',E13.6,10X,'RHSW7 = ',E13.6)
!-----------------------------------------------------------------------
!     End of subroutine PATH
!-----------------------------------------------------------------------
      END
      SUBROUTINE DRAIN2P(HAW,PHAW,I)
!-----------------------------------------------------------------------
!     Called by PATH
!
!     This subroutine calculates the apparent water saturation, 
!     entrapped air saturation, permeability, and apparent saturation
!     capacity for the primary drainage path with and without air
!     entrapment. 
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
!
!     Effective apparent water saturation
!
      ASW = (1.D0+(ALPHAD*HAW)**XN)**(-XM)
      IF (ASW .GT. 1.0) ASW = 1.0
      IF (ASW .LT. RSW) THEN
        IF (HAW .GT. HYSMXH) THEN
          RSW = ASWMXH(I)
        ELSE
          RSW = ASW
        ENDIF
        IF (RAW .GT. 0.0) SRAW = (1.D0-RSW)/(1.D0+RAW*(1.D0-RSW))
        IF (SRAW .LT. AIRTOL) SRAW = 0.0
      ENDIF
      SATW = 0.0
      IF (HAW .LE. HYSMXH) THEN
        IF (RSW .LT. 1.0) SATW = SRAW*(ASW-RSW)/(1.-RSW)
      ENDIF
!
!     Effective and actual saturations for all options.
!
      ESW  = ASW-SATW
      SW   = ESW*(1.D0-SM)+SM
      ASAT = SATW*(1.D0-SM)
      X1   = ASW**XXM
!
!     Calculate capacity (CWW)
!
      CWW  = PHI*(1.D0-SM)*ALPHAD*(XN-1.D0)*X1*(1.D0-X1)**XM
!
!     Note: Lenhard used ESW for PO calculation.  This is not consistent
!     with his using ASW in DRY2P and WET2P.  Using ASW increases K.
!
      P0 = 1.D0-(1.D0-ESW**XXM)**XM
      P3 = 0.0
      IF (IPATH .EQ. 1) THEN
!
!     Corrections to permeability and capacity for entrapped air 
!     effects in the absence of hysteresis 
!
        IF (SATW .GT. 0) THEN
          CSATW = SRAW/(1.D0-RSW)
          X7  = 1.0D0-CSATW
          CWW = X7*CWW
          P3  = CSATW*((1.D0-RSW**XXM)**XM-(1.D0-ASW**XXM)**XM)
        ENDIF
      ENDIF
!
!     Calculate relative permeability (PERMW)
!
      PERMW = (ESW**0.5)*(P0-P3)**2
!
!     Check for unusual permeability or capacity values
!
      IF (PERMW .LT. 0.0 .OR. PERMW .GT. 1.0) WRITE(LUS,3) I,PERMW
      IF (CWW .LT. 0.0) WRITE(LUS,4) I,CWW
      RETURN
!-----------------------------------------------------------------------
!     FORMAT Statements
!-----------------------------------------------------------------------
3     FORMAT(' Relative permeability in DRAIN2P is < 0.0 or > 1.0',/,   &
     &' NODE = ',I5,', PERMW = ',1PE15.7)
4     FORMAT(' Water capacity in DRAIN2P is negative',/,                &
     &' NODE = ',I5,', CWW = ',1PE15.7)
!-----------------------------------------------------------------------
!     End of subroutine DRAIN2P
!-----------------------------------------------------------------------
      END
      SUBROUTINE DRY2P (HAW,PHAW,I)
!-----------------------------------------------------------------------
!     Called by PATH
!
!     This subroutine calculates water and air saturations,
!     relative permeability, and capacity for drainage paths
!
!     Variable list:
!
!     RDISWI = Reversal point from the previous, outer drainage path
!              to the imbibition curve connecting it to the current
!              inner drainage path
!     RIDSWI = Reversal point on the most recent imbibition curve to the
!              current drainage path
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
!
!     Calculate effective apparent water saturation
!
      RIDSWD = (1.D0+(ALPHAD*RHSW(I,IPSW1))**XN)**(-XM)
      RDISWD = (1.D0+(ALPHAD*RHSW(I,IPSW1-1))**XN)**(-XM)
      RIDRDI = RIDSWD-RDISWD
      SWD    = (1.D0+(ALPHAD*HAW)**XN)**(-XM)
      X8     = 1.0
      IF (RIDRDI .GT. 0.0) X8 = (RASW(I,IPSW1)-RASW(I,IPSW1-1))/RIDRDI
      IF (RIDRDI .LT. 0.0) THEN
        WRITE(LUS,*) ' STOP: RIDRDI calculated to be < 0.0 in DRY2P'
        WRITE(LUS,5) I,IPSW1,JJJH,HAW,SWD,RHSW(I,IPSW1-1),RHSW(I,IPSW1),&
     &    RASW(I,IPSW1-1),RASW(I,IPSW1),RIDSWD,RDISWD
        GO TO 999
      ENDIF
      ASW = X8*(SWD-RDISWD)+RASW(I,IPSW1-1)
      IF (ASW .GT. 1.0) ASW = 1.0
      SATW = 0.0
      IF (RSW .LT. 1.0) SATW = SRAW*(ASW-RSW)/(1.D0-RSW)
      ESW  = ASW-SATW
      SW   = ESW*(1.D0-SM)+SM
      ASAT = SATW*(1.D0-SM)
      P0   = 1.D0-((1.D0-ASW**XXM)**XM)
!
!     Calculate corrections to relative permeability and capacity from
!     entrapped air effects
!
      P3 = 0.0
      X7 = 1.0
      CSATW = 0.0
      IF (SATW .GT. 0.0) THEN
        CSATW = SRAW/(1.D0-RSW)
        X7 = 1.D0-CSATW
        P3 = CSATW*((1.D0-RSW**XXM)**XM-(1.D0-ASW**XXM)**XM)
      ENDIF
!
!     Calculate relative permeability
!
      PERMW = (ESW**.5)*(P0-P3)**2.
      X1 = SWD**XXM
!
!     Calculate capacity
!
      CWW = PHI*(1.D0-SM)*ALPHAD*(XN-1.D0)*X1*X7*X8*(1.D0-X1)**XM
      IF (PERMW .LT. 0.0 .OR. PERMW .GT. 1.0) WRITE(LUS,3) I,PERMW
      IF (CWW .LT. 0.0) WRITE(LUS,4) I,CWW
      RETURN
999   STOP
!-----------------------------------------------------------------------
!     FORMAT Statements
!-----------------------------------------------------------------------
3     FORMAT(' Relative permeability in DRY2P is < 0.0 or > 1.0',/,     &
     &' NODE = ',I5,', PERMW = ',1PE15.7)
4     FORMAT(' Water capacity in DRY2P is negative',/,                  &
     &' NODE = ',I5,', CWW = ',1PE15.7)
5     FORMAT(/,' I = ',I2,'  IPSW1 = ',I2,'  JJJH = ',I2,/,1P,          &
     &' HAW             = ',E15.7,'  SWD           = ',E15.7,/,         &
     &' RHSW(I,IPSW1-1) = ',E15.7,'  RHSW(I,IPSW1) = ',E15.7,/,         &
     &' RASW(I,IPSW1-1) = ',E15.7,'  RASW(I,IPSW1) = ',E15.7,/,         &
     &' RIDSWD          = ',E15.7,'  RDISWD        = ',E15.7)
!-----------------------------------------------------------------------
!     End of subroutine DRY2P
!-----------------------------------------------------------------------
      END
      SUBROUTINE WET2P (HAW,PHAW,I)
!-----------------------------------------------------------------------
!     Called by PATH
!
!     This subroutine calculates water and entrapped air saturations,
!     relative permeability, and capacity for wetting paths.
!
!     RDISWI = Reversal point on the saturation imbibition curve from
!              the primary drainage path to the imbibition path
!     RIDSWI = Reversal point on the saturation imbibition curve from
!              the imbibition path to the primary drainage path
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
!
!     Calculate effective apparent water saturation
!
      RIDSWI = (1.D0+(ALPHAI*RHSW(I,IPSW1-1))**XN)**(-XM)
      RDISWI = (1.D0+(ALPHAI*RHSW(I,IPSW1))**XN)**(-XM)
      RIDRDI = RIDSWI-RDISWI
      SWI    = (1.D0+(ALPHAI*HAW)**XN)**(-XM)
      IF (RIDRDI .GT. 0) THEN
        X8 = (RASW(I,IPSW1-1)-RASW(I,IPSW1))/RIDRDI
      ELSE IF (RIDRDI .EQ. 0) THEN
        X8 = 1.0
      ELSE
        WRITE(LUS,*) ' STOP: RIDRDI calculated to be < 0.0 in WET2P'
        WRITE(LUS,5) I,IPSW1,JJJH,HAW,SWI,RHSW(I,IPSW1-1),RHSW(I,IPSW1),&
     &    HAW,SWI,RASW(I,IPSW1-1),RASW(I,IPSW1),RIDSWI,RDISWI
        GO TO 999
      ENDIF
      ASW = X8*(SWI-RDISWI)+RASW(I,IPSW1)
      IF (ASW .GT. 1.0) ASW = 1.0
      SATW = 0.0
      IF (RSW .LT. 1.0) SATW = SRAW*(ASW-RSW)/(1.D0-RSW)
      ESW  = ASW-SATW
      SW   = ESW*(1.D0-SM)+SM
      ASAT = SATW*(1.D0-SM)
      P0 = 1.D0-((1.D0-ASW**XXM)**XM)
      P3 = 0.0
      X7 = 1.0
      IF (SATW .GT. 0.0) THEN
!
!     Calculate correction terms for air entrapment effects on capacity
!     and relative permeability
!
        CSATW = SRAW/(1.D0-RSW)
        X7    = 1.D0-CSATW
        P3    = CSATW*((1.D0-RSW**XXM)**XM-(1.D0-ASW**XXM)**XM)
      ENDIF
!
!     Calculate relative permeability
!
      PERMW = (ESW**.5)*(P0-P3)**2.
      X1 = SWI**XXM
!
!     Calculate capacity
!
      CWW = PHI*(1.D0-SM)*ALPHAI*(XN-1.D0)*X1*X8*X7*(1.D0-X1)**XM
      IF (PERMW .LT. 0.0 .OR. PERMW .GT. 1.0) WRITE(LUS,3) I,PERMW
      IF (CWW .LT. 0.0) WRITE(LUS,4) I,CWW
      RETURN
999   STOP
!-----------------------------------------------------------------------
!     FORMAT Statements
!-----------------------------------------------------------------------
3     FORMAT(' Relative permeability in WET2P is < 0.0 or > 1.0',/,     &
     &' NODE = ',I5,', PERMW = ',1PE15.7)
4     FORMAT(' Water capacity in WET2P is negative',/,                  &
     &' NODE = ',I5,', CWW = ',1PE15.7)
5     FORMAT(/,' I = ',I2,'  IPSW1 = ',I2,'  JJJH = ',I2,/,1P,          &
     &' HAW             = ',E15.7,'            SWI = ',E15.7,/,         &
     &' RHSW(I,IPSW1-1) = ',E15.7,'  RHSW(I,IPSW1) = ',E15.7,/,         &
     &' RASW(I,IPSW1-1) = ',E15.7,'  RASW(I,IPSW1) = ',E15.7,/,         &
     &' RIDSWI          = ',E15.7,'  RDISWI        = ',E15.7)
!-----------------------------------------------------------------------
!     End of subroutine WET2P
!-----------------------------------------------------------------------
      END
      SUBROUTINE UPDATE(HAW,PHAW,I)
!-----------------------------------------------------------------------
!     Called by UPDATE
!
!     Reset global variables after successful solution.
!     Check the index variables and set them to the proper saturation
!     path.
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (ASW .LT. RSWAW(I)) RSWAW(I) = ASW
      IF (HAW .GT. HYSMXH) RSWAW(I) = ASWMXH(I)
      SARW(I) = SRAW
      IF (IPATH .EQ. 1) RETURN
!
!     Hysteresis is operative
!
      IF (HAW .EQ. 0.0 .AND. IPSW1 .GT. 1) THEN
!
!       node is saturated; reset IPSW, JJH, and RASW/RHSW (I,3:IPATH)
!      
        IPSW1 = 2
        JJJH  = 1
        RASW(I,3) = 1.0
        RHSW(I,3) = 0.0
        DO M = 5,IPATH,2
          RASW(I,M) = ASWMXH(I)
          RHSW(I,M) = HYSMXH
        ENDDO
        DO M = 4,IPATH,2
          RASW(I,M) = 1.0
          RHSW(I,M) = 0.0
        ENDDO
      ENDIF
!
!       escaped from an inner hysteretic path to an outer path;
!       reset RASW/RHSW values for paths higher than IPSW1
!     
      IF (LSW .GT. 0) THEN
        DO M = IPSW1+1,IPATH
          IF ((M-(M/2)*2) .EQ. 0) THEN
!
!         drying to wetting path reversal point
!
            RASW(I,M) = 1.0
            RHSW(I,M) = 0.0
          ELSE
!
!         wetting to drying path reversal point
!
            RASW(I,M) = ASWMXH(I)
            RHSW(I,M) = HYSMXH
          ENDIF
        ENDDO
      ENDIF
!
!     reset variables to primary drainage path if HAW >= RHSW(I,2)
!
      IF (HAW .GE. RHSW(I,2)) THEN
        IPSW1 = 1
        JJJH  = 0
        IF (HAW .GT. HYSMXH) THEN
          RASW(I,2) = ASWMXH(I)
          RHSW(I,2) = HYSMXH
        ELSE
          RASW(I,2) = ASW
          RHSW(I,2) = HAW
        ENDIF
        DO M = 3,IPATH,2
          RASW(I,M) = ASWMXH(I)
          RHSW(I,M) = HYSMXH
        ENDDO
        DO M = 4,IPATH,2
          RASW(I,M) = 1.0
          RHSW(I,M) = 0.0
        ENDDO
      ENDIF
      JJH(I) = JJJH
      IPSW(I) = IPSW1
      IF (IPSW1 .EQ. IPATH) RETURN
!
!     check reversal points and reset to current values of HAW and ASW
!
      IF ((JJJH .EQ. 0 .AND. HAW .GT. RHSW(I,IPSW1+1)) .OR.             &
     &    (JJJH .EQ. 1 .AND. HAW .LT. RHSW(I,IPSW1+1))) THEN
        RHSW(I,IPSW1+1) = HAW
        RASW(I,IPSW1+1) = ASW
        IF (HAW .GT. HYSMXH) THEN
          RHSW(I,IPSW1+1) = HYSMXH
          RASW(I,IPSW1+1) = ASWMXH(I)
        ENDIF
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine UPDATE
!-----------------------------------------------------------------------
      END
      SUBROUTINE HYSOUT
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     This subroutine outputs hysteresis variables from UNSAT-H to
!     an external file with the input filename, and extension *.HRI,
!     where HRI stands for Hysteresis Restart Information.
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      CHARACTER*80 HFILE
      HFILE = RFILE
      NCHRH = INDEX(HFILE,'.')
      HFILE(NCHRH:NCHRH+3) = '.hri'
      OPEN(UNIT=LUH,FILE=HFILE,STATUS='UNKNOWN')
      CLOSE(UNIT=LUH,STATUS='DELETE')
      OPEN(UNIT=LUH,FILE=HFILE,STATUS='NEW')
!
!     Write initial and reversal values to restart file
!
      WRITE(LUH,2005) (IPSW(I),I=1,NPT)
      WRITE(LUH,2005)  (JJH(I),I=1,NPT)
      WRITE(LUH,2010) (SARW(I),I=1,NPT)
      WRITE(LUH,2010)(RSWAW(I),I=1,NPT)
      DO I=1,NPT
        IPATH = IPATHA(MAT(I))
        WRITE(LUH,2010)(RASW(I,J),J=1,IPATH)
        WRITE(LUH,2010)(RHSW(I,J),J=1,IPATH)
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     FORMAT Statements
!-----------------------------------------------------------------------
2005  FORMAT(1X,70I1)
2010  FORMAT(1P7E15.7)
!-----------------------------------------------------------------------
!     End of subroutine HYSOUT
!-----------------------------------------------------------------------
      END
      SUBROUTINE LOWBC(IPRINT,LUSLUW)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     Read in the lower boundary conditions
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION NQLEAK(366),QLEAK(366)
      IF (IPRINT .GT. 0) WRITE(LOUT,'(80("-"),/,A)')                    &
     &' Lower Boundary Option:'
      IF (LOWER .EQ. 1) THEN
        IF (IPRINT .GT. 0) WRITE(LOUT,7010) LOWER
      ELSE IF (LOWER .EQ. 2) THEN
        IF (IPRINT .GT. 0) WRITE(LOUT,7020) LOWER
      ELSE IF (LOWER .EQ. 3) THEN
        IF (IPRINT .GT. 0) WRITE(LOUT,7030) LOWER
        READ(LUI,'(I)') NTLEAK
        READ(LUI,'(5(I,D))') (NQLEAK(I),QLEAK(I),I=1,NTLEAK)
        IF (NTLEAK .GT. 366) THEN
          WRITE(LUS,'(A)') 'NQLEAK,QLEAK dimensions exceed 366'
          GOTO 999
        ENDIF
        IF (IPRINT .GT. 2) THEN
          WRITE(LOUT,7310) NTLEAK
          WRITE(LOUT,7320) (NQLEAK(I),QLEAK(I),I=1,NTLEAK)
        ENDIF
        JJ = 1
        DO I=1,366
          IF (I .EQ. NQLEAK(JJ)) THEN
            QWLEAK(I) = QLEAK(JJ)
            JJ = JJ+1
            J  = JJ-1
          ELSE
            IF (I .GT. NQLEAK(NTLEAK)) THEN
              QWLEAK(I) = 0.0
            ELSE
              QWLEAK(I) = QLEAK(J)+(I-NQLEAK(J))*(QLEAK(JJ)-QLEAK(J))/  &
     &        (NQLEAK(JJ)-NQLEAK(J))
            ENDIF
          ENDIF
        ENDDO
        IF (IPRINT .GT. 2) THEN
          WRITE(LOUT,7330)
          WRITE(LOUT,7320) (I,QWLEAK(I),I=1,366)
        ENDIF
      ELSE IF (LOWER .EQ. 4) THEN
        IF (IPRINT .GT. 0) WRITE(LOUT,7040) LOWER
      ELSE IF (LOWER .EQ. 5) THEN
        IF (IPRINT .GT. 0) WRITE(LOUT,7050) LOWER
        READ(LUI,'(I)') INDEXB
        DO I=1,INDEXB
          READ(LUI,'(2D)') BBHEAD(I,1),BBHEAD(I,2)
        ENDDO
        IF (IPRINT .GT. 2) THEN
          WRITE(LOUT,7335) INDEXB
          WRITE(LOUT,7336) ((I,BBHEAD(I,1),BBHEAD(I,2)),I=1,INDEXB)
        ENDIF
      ELSE
        WRITE(LUS,7340) LOWER
        GOTO 999
      ENDIF
      RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in LOWBC'
      IF (LUSLUW) WRITE(LUS,'(A)')                                      &
     &' Program terminated prematurely in LOWBC'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
7010  FORMAT(/,T4,'LOWER =',I2,':  unit gradient')
7020  FORMAT(/,T4,'LOWER =',I2,':  constant head')
7030  FORMAT(/,T4,'LOWER =',I2,':  user-supplied flux variations')
7040  FORMAT(/,T4,'LOWER =',I2,':  impermeable')
7050  FORMAT(/,T4,'LOWER =',I2,':  user-supplied head variations')
7300  FORMAT(I5/,100(5(I5,F5.0)/))
7310  FORMAT(/,T4,'NTLEAK =',I3,' flux changes',//,                     &
     &T4'Input Values (cm/d)')
7320  FORMAT(T4,5(' Day   Flux '),/,T4,5(' ---   ---- '),/,             &
     &(T4,5(I4,F8.4)))
7330  FORMAT(/,T4,'Output values (cm/d) interpolated from input values')
7335  FORMAT(/,T4'INDEXB =',I3,' pairs of time-head values',/,          &
     &T4'Index',T14,'Time (d)',T29,'Head (cm)')
7336  FORMAT(I7,1P2E15.7)
7340  FORMAT(/,' ERROR:  LOWER options are 1 to 5; input value was',I2)
!-----------------------------------------------------------------------
!     End of subroutine LOWBC
!-----------------------------------------------------------------------
      END
      SUBROUTINE PLANTIN(IPRINT,LUSLUW,MATN)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     Reads in plant parameters
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION IDLAI(366),VLAI(366)
      IF(IPRINT .GT. 0) WRITE(LOUT,'(80("-"),/,A)') ' Plant parameters:'
!
!     Read LEAF,NFROOT,NUPTAK,NFPET,NSOW,NHRVST
!
      READ(LUI,'(6I)') LEAF,NFROOT,NUPTAK,NFPET,NSOW,NHRVST
      READ(LUI,'(D)') BARE
      IF (IPRINT .GT. 1) WRITE(LOUT,6000) LEAF,NFROOT,NUPTAK,NFPET,     &
     &    NSOW,NHRVST,BARE
!
!     LEAF = 1   Leaf area index IDLAI,VALLAI
!          = 2   User-supplied subroutine for leaf area index
!
      IF (LEAF .EQ. 1) THEN
        READ(LUI,'(I)') NDLAI
        IF (NDLAI .GT. 366) THEN
          WRITE(LOUT,'(A)') ' NDLAI > 366, VLAI dimension exceeded'
          GOTO 999
        ENDIF
        READ(LUI,'(4(I,D))') (IDLAI(I),VLAI(I),I=1,NDLAI)
        IF (IPRINT .GT. 2) THEN
          WRITE(LOUT,6020) NDLAI,(IDLAI(I),VLAI(I),I=1,NDLAI)
          WRITE(LOUT,6025)
        ENDIF
        J = 1
        DO I=1,NDAYS
          IF (I .LT. IDLAI(1)) THEN
            DLAI(I) = 0.0
          ELSE IF (I .GT. IDLAI(NDLAI)) THEN
            DLAI(I)=0.0
          ELSE IF(I .GE. IDLAI(1) .AND. I .LE. IDLAI(NDLAI)) THEN
            IF (I .GE. IDLAI(J) .AND. J .LT. NDLAI) J = J+1
            JJ      = J-1
            DIFF    = (VLAI(J)-VLAI(JJ))/(IDLAI(J)-IDLAI(JJ))
            DLAI(I) = VLAI(JJ)+DIFF*(I-IDLAI(JJ))
          ENDIF
        ENDDO
        IF(IPRINT.GT.2) WRITE(LOUT,6030) (I,DLAI(I),I=1,NDAYS)
      ELSE IF (LEAF .EQ. 2) THEN
        WRITE(LUS,'(A)') ' Subroutine MYLAI was not provided'
        GOTO 999
      ELSE IF (LEAF .GT. 2) THEN
        WRITE(LUS,6040) LEAF
        IF (LUSLUW) WRITE(LOUT,6040) LEAF
        GOTO 999
      ENDIF
      IF (NFROOT .NE. 1 .AND. NFROOT .NE. 2) THEN
        WRITE(LUS,6050) NFROOT
        IF (LUSLUW) WRITE(LOUT,6050) NFROOT
        GOTO 999
      ENDIF
      IF (NFROOT .EQ. 1) THEN
!
!     NFROOT = 1:  Read AA,B1,B2,NTROOT
!
        READ (LUI,'(3D)') AA,B1,B2
        IF (IPRINT .GT. 1) WRITE(LOUT,6060) AA,B1,B2
        READ (LUI,'(10I)') (NTROOT(I),I=1,NPT)
         NGDAYS = NHRVST-NSOW+1
        IF (NSOW .GT. NHRVST) NGDAYS = ILEAP+365-NSOW+1+NHRVST
        I = 1
630     IF (NTROOT(I) .LT. NGDAYS .AND. I .LE. NPT) THEN
          I = I+1
          GO TO 630
        ENDIF
        MDEPTH = I-1
        MXROOT = MDEPTH
        IF (IPRINT .GT. 1) THEN
          WRITE(LOUT,6070)
          WRITE(LOUT,6080) 
!
!     Estimate root density factor (RDF) at each node of root depth
!
          ROOTS  = 0.0
          ZZ0    = 0.0
          RDF(1) = 0.0
          DO I=2,MDEPTH
            RDF(I) = RLD(Z(I))
            ROOTS  = ROOTS+RDF(I)*(ZZ(I)+ZZ0)*0.5
            ZZ0    = ZZ(I)
          ENDDO
          ROOTS = ROOTS-RDF(MDEPTH)*ZZ0*0.5
          WRITE(LOUT,6090)                                              &
     &      (NTROOT(I),Z(I),RDF(I),RDF(I)/ROOTS,I=1,MDEPTH)
        ENDIF
        IF (IPRINT .GT. 0) THEN
          WRITE(LOUT,6100) MXROOT
        ENDIF
!
!     NFROOT = 2   User-defined root density
!
      ELSE
        WRITE(LUS,'(A)') ' Subroutine "MYROOT" was not provided'
        WRITE(LOUT,'(A)') ' Subroutine "MYROOT" was not provided'
      ENDIF
      IF (NUPTAK .EQ. 1) THEN
!
!     Read HW,HD,HN if NUPTAK = 1
!
        IF (IPRINT .GT. 1) WRITE(LOUT,6110)
        MATSAV = MAT(1)
        DO I=1,MATN
          MAT(1) = I
          READ (LUI,'(3D)') HW,HD,HN
          IF (HW .LE. HD) THEN
            WRITE(LOUT,'(A,1P2E15.7)') ' HW and HD = ',HW,HD
            WRITE(LOUT,*) ' Error:  HW must be greater than HD'
            GOTO 999
          ENDIF
          IF (HD .LE. HN) THEN
            WRITE(LOUT,'(A,1P2E15.7)') ' HD and HN = ',HD,HN
            WRITE(LOUT,*) ' Error:  HD must be greater than HN'
            GOTO 999
          ENDIF
          CALL RETENT(1,1,HW,1,THETAW(I),0,C)
          CALL RETENT(1,1,HD,1,THETAD(I),0,C)
          CALL RETENT(1,1,HN,1,THETAN(I),0,C)
          IF (IPRINT .GT. 1) WRITE(LOUT,6120)I,THETAW(I),THETAD(I),     &
     &        THETAN(I)
        ENDDO
        MAT(1) = MATSAV
      ELSE IF (NUPTAK .EQ. 2) THEN
        WRITE(LUS,'(A)') ' User subroutine "MYSINK" was not provided'
        IF (LUSLUW)                                                     &
     &    WRITE(LOUT,'(A)') ' User subroutine "MYSINK" was not provided'
        GOTO 999
      ELSE
        WRITE(LUS,6130) NUPTAK
        IF (LUSLUW) WRITE(LOUT,6130) NUPTAK
        GOTO 999
      ENDIF
      IF (NFPET .EQ. 1) THEN
        READ (LUI,'(5D)') (PETPC(I),I=1,5)
      ELSEIF (NFPET .EQ. 2) THEN
        READ (LUI,'(D)') BIOMAS
      ENDIF
      RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in PLANTIN'
      IF (LUSLUW)                                                       &
     &  WRITE(LOUT,'(A)') ' Program terminated prematurely in PLANTIN'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
6000  FORMAT(/,'   LEAF=',I3,', NFROOT=',I3,', NUPTAK=',I3,             &
     &', NFPET=',I3,', NSOW=',I3,', NHRVST=',I3,', BARE=',F5.3)
6020  FORMAT(/,'   Number of Growth Day (NDLAI) - Leaf Area Index ',    &
     &'(VLAI) pairs',I5,//,T12,'NDLAI',T22,'VLAI',/,T12,5('-'),         &
     &T22,5('-'),/,(T12,I5,T21,F6.3))
6025  FORMAT(/,6('   DAY   LAI '),/,6('   ---  -----'))
6030  FORMAT(6(I6,F7.3))
6040  FORMAT(/,' ERROR:  LEAF options are 0 and 1. Input value was',I4)
6050  FORMAT(/,' ERROR:  NFROOT options are 1. Input value was ',I4)
6060  FORMAT(80('-'),/,' NFROOT = 1:',                                  &
     &'   Negative exponential representation of root growth',/,        &
     &'   AA (intersection of curve at z=0 with abscissa)    =',F10.3,/,&
     &'   B1 (coefficient defining degree of curvature)      =',F10.3,/,&
     &'   B2 (coefficient determining the value of asymptote =',F10.3)
6070  FORMAT(/,'   Root depth, density, and weight/node versus depth',/)
6080  FORMAT(T4,'DAY',T14,'MAX',T25,'ROOT',T33,'NORMALIZED',/,          &
     &    T10,'ROOT DEPTH',T23,'DENSITY',T35,'DENSITY',/,               &
     &    T23,'(cm/cm)',T35,'(1/cm)',/,                                 &
     &    T4,'---',T10,'----------',T23,'-------',T33,'----------')
6090  FORMAT(T4,I3,T12,F7.2,T20,F9.3,T34,F7.4)
6100  FORMAT(/,                                                         &
     &'   MXROOT (deepest node that roots penetrate) =',I4,/,80("-"))
6110  FORMAT(' NUPTAK = 1:',                                            &
     &'  Feddes et al. 1975 moisture dependent sink term')
6120  FORMAT(/'   For Material No.',I4,/,                               &
     &'    THETAW (wilting point moisture content)          = ',F10.4,/,&
     &'    THETAD (lower limit of optimum moisture content) = ',F10.4,/,&
     &'    THETAN (upper limit of optimum moisture content) = ',F10.4)  &
6130  FORMAT(/' ERROR:  NUPTAK options are 1. Input value was ',I4)
!-----------------------------------------------------------------------
!     End of subroutine PLANTIN
!-----------------------------------------------------------------------
      END
      SUBROUTINE ETPAR(IPRINT,LUSLUW,NFHOUR)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     Calls MYHRLY
!
!     Reads in PET/MET parameters
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (IPRINT .GT. 0) THEN
        WRITE(LOUT,'(80("-"),/,A,/)') ' ET parameters:'
        IF (IEVOPT .EQ. 0 .AND. IPLANT .EQ. 0) THEN
          WRITE(LOUT,'(A)') '   IEVOPT = 0 and IPLANT = 0 (no ET)'
          WRITE(LOUT,'(A)') '   Thus, no ET parameters required'
          GOTO 998
        ENDIF
      ENDIF
!
!     Read FPET (NFHOUR = 1) or calculate FPET (NFHOUR = 2)
!
      IF ((IEVOPT .EQ. 1 .OR. IPLANT .EQ. 1) .AND. IHEAT .EQ. 0) THEN
        IF (NFHOUR .EQ. 1) THEN
          READ(LUI,'(8D)',ERR=999) (FPET(I),I=1,24)
          IF (IPRINT .GT. 1) WRITE(LOUT,7120) 
        ELSE IF (NFHOUR .EQ. 2) THEN
          CALL MYHRLY(FPET)
          IF (IPRINT .GT. 1) WRITE(LOUT,7130)
        ELSE
          WRITE(LOUT,7140) NFHOUR
          GOTO 999
        ENDIF
        SUM = 0.0
        DO I=1,24
          SUM = SUM+FPET(I)
        ENDDO
        IF (SUM .LT. 0.999 .OR. SUM .GT. 1.001) THEN
          WRITE(LOUT,7150) SUM
          WRITE(LOUT,'(3(/8(3X,F6.4)))') (FPET(I),I=1,24)
          GOTO 999 
        ELSE
          DO I=1,24
            FPET(I) = FPET(I)/SUM
          ENDDO
          IF (IPRINT .GT. 1) THEN
            WRITE(LOUT,'(3(/8(3X,F6.4)))') (FPET(I),I=1,24)
          ENDIF
        ENDIF
      ENDIF
7120  FORMAT('   NFHOUR = 1:  Hourly PET distribution values provided')
7130  FORMAT('   NFHOUR = 2:  User subroutine for hourly PET provided')
7140  FORMAT(/,' ERROR:  Current NFHOUR options are 1 and 2.',/,        &
     &' Input value was ',I2) 
7150  FORMAT(' ERROR:  Sum of diurnal PET distribution factors is ',    &
     &1PE15.7,'.  It should be 1.0')
!
!     Read MET parameters for PET calc.(IHEAT=0) or heat flow (IHEAT=1)
!
      IF (IETOPT .EQ. 1) THEN
        IF (IHEAT .EQ. 0) THEN
          READ(LUI,'(4D)',ERR=999) ALBEDO,ALT,ZU,PMB
          IF (IPRINT .GT. 1) WRITE(LOUT,7166) ALBEDO,ALT,ZU,PMB
        ELSE
          IF (UPPERH .EQ. 0) THEN
            READ(LUI,'(6D)',ERR=999) ZH,ZM,ZT,ZU,D,LAT
            IF (IPRINT .GT. 1) WRITE(LOUT,7167) ZH,ZM,ZT,ZU,D,LAT
            DAYSEXT = SEXT*86400.0
            TPIY   = 2.0*PI/365.0
            LAT    = LAT*2.0*PI/360.0
            SINLAT = SIN(LAT)
            COSLAT = COS(LAT)
            TANLAT = TAN(LAT)
          ENDIF
        ENDIF
      ENDIF
7166  FORMAT(' IETOPT = 1 and IHEAT = 0:',/,
     &' PET calculated from meteorological data',/,                     &
     &' using subroutine CALPET',//,1P,                                 &
     &' ALBEDO = ',E10.3,/,                                             &
     &'    ALT = ',E10.3,' (m)',/,                                      &
     &'     ZU = ',E10.3,' (m)',/,                                      &
     &'    PMB = ',E10.3,' (mb)',/)
7167  FORMAT(' IHEAT = 1:',/,                                           &
     &' Evaporation calculated from meterological data and heat',/,     &
     &' flow calculations',/,1P,                                        &
     &'  ZH = ',E10.3,' (m)',/,                                         &
     &'  ZM = ',E10.3,' (m)',/,                                         &
     &'  ZT = ',E10.3,' (m)',/,                                         &
     &'  ZU = ',E10.3,' (m)',/,                                         &
     &'   D = ',E10.3,' (m)',/,                                         &
     &' LAT = ',E10.3,' (deg)',/)
998   RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in ETPAR'
      IF (LUSLUW) WRITE(LOUT,'(A)')                                     &
     &' Program terminated prematurely in ETPAR'
      STOP
!-----------------------------------------------------------------------
!     End of subroutine ETPAR
!-----------------------------------------------------------------------
      END
      SUBROUTINE ETBC(IPRINT,LUSLUW,LULI)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     Reads in the PET/MET boundary conditions
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (IETOPT .EQ. 0) THEN
        READ(LULI,'(8D)',ERR=999) (PET(I),I=1,NDAYS)
      ELSE 
        READ(LULI,*,ERR=999) ((RMDATA(I,J),I=1,8),J=1,NDAYS)
        IF (IPRINT .GT. 2) WRITE(LOUT,7168) ((RMDATA(I,J),              &
     &      I=1,8),J=1,NDAYS)
        IF (IHEAT .EQ. 0)  THEN
          CALL CALPET
!
!     CALPET utilizes data with same units as FAOPET
!
!     Once through CALPET, data are converted to UNSAT-H units
!     RMDATA(2) = TMAX (F) converted to K
!     RMDATA(3) = TMIN (F) converted to K
!     RMDATA(4) = TDEW (F) converted to K
!     RMDATA(8) = precipitation (in.) converted to cm
!
          DO J=1,NDAYS
            RMDATA(2,J) = 273.16+(RMDATA(2,J)-32.0)*5.0/9.0
            RMDATA(3,J) = 273.16+(RMDATA(3,J)-32.0)*5.0/9.0
            RMDATA(4,J) = 273.16+(RMDATA(4,J)-32.0)*5.0/9.0
            RMDATA(8,J) = RMDATA(8,J)*2.54
          ENDDO
        ELSE
!-----------------------------------------------------------------------
!     RMDATA(2) = TMAX (F) converted to K
!     RMDATA(3) = TMIN (F) converted to K
!     RMDATA(4) = TDEW (F) converted to K
!     RMDATA(5) = SR_MEAS (ly/d) converted to J/m2/d
!     RMDATA(6) = WIND (mph) converted to m/s
!     RMDATA(7) = CLOUD (tenths) converted to fraction
!     RMDATA(8) = precipitation (in.) converted to cm
!-----------------------------------------------------------------------
          DO J=1,NDAYS
            RMDATA(2,J) = 273.16+(RMDATA(2,J)-32.0)*5.0/9.0
            RMDATA(3,J) = 273.16+(RMDATA(3,J)-32.0)*5.0/9.0
            RMDATA(4,J) = 273.16+(RMDATA(4,J)-32.0)*5.0/9.0
            RMDATA(5,J) = RMDATA(5,J)*41880.0
            RMDATA(6,J) = RMDATA(6,J)*0.447
            RMDATA(7,J) = RMDATA(7,J)/10.0
            RMDATA(8,J) = RMDATA(8,J)*2.54
          ENDDO
        ENDIF
      ENDIF
      IF (IHEAT .EQ. 0) THEN
        IF (IPLANT .EQ. 0 .AND. IEVOPT .GT. 0) THEN
          DO I=1,NDAYS
            PEVAPO(I) = PET(I)
            TOTPEV = TOTPEV+PEVAPO(I)
          ENDDO
          IF (IPRINT .GT. 0) THEN
            WRITE(LOUT,'(80("-"),//,A)')                                &
     &      ' Potential Evaporation Values (no plants):'
            IF (IPRINT .GT. 1) THEN
              WRITE(LOUT,7170)
              WRITE(LOUT,7180) (I,PEVAPO(I),I=1,NDAYS)
            ENDIF
            WRITE(LOUT,7190) TOTPEV
          ENDIF
        ENDIF
      ENDIF
      RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in ETBC'
      IF (LUSLUW) WRITE(LOUT,'(A)')                                     &
     &' Program terminated prematurely in ETBC'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
7168  FORMAT(/,' IETOPT = 1:  Meteorological Data',//,                  &
     &    '          Temperature   Solar   Wind  Cloud',/,              &
     &    '  IDAY  Max  Min   Dew   Rad   Speed  Cover  Prec.'/,        &
     &    '         F    F     F    ly/d   mph   tenth    in '/,        &
     &    '  ----  ---  ---  ----  -----  -----  -----  -----'/,        &
     &    (F6.0,F5.0,F5.0,F6.1,F7.0,2F7.1,F7.2))
7170  FORMAT(/,(5(' DAY  PEVAPO   '),/,5(' ---  ------   ')))
7180  FORMAT(5(I4,F8.4,3X))
7190  FORMAT(/' Totals:  PEVAPO = ',F8.4)
!-----------------------------------------------------------------------
!     End of subroutine ETBC
!-----------------------------------------------------------------------
      END
      SUBROUTINE PETPART(IPRINT,LUSLUW)
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     Partitions PET into PEVAPO and PTRANS
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF(IPRINT .GT. 0)                                                 &
     &WRITE(LOUT,'(80("-"),/,A)') ' PET partitioning:'
      IF (NFPET .EQ. 1) THEN
!
!     If NFPET = 1, Ritchie method of partitioning PET 
!
        PCA = PETPC(1)
        PCB = PETPC(2)
        PCC = PETPC(3)
        PCD = PETPC(4)
        PCE = PETPC(5)
        IF (IPRINT .GT. 1) THEN
          WRITE(LOUT,7200) PCA,PCB,PCC,PCD,PCE
          WRITE(LOUT,7210)
        ENDIF
        DO LDAY=1,NDAYS
          IF (NSOW .LE. LDAY) THEN
            NGROWD = LDAY-NSOW+1
          ELSE
            NGROWD = ILEAP+365+LDAY-NSOW+1
          ENDIF
          IF (NGROWD .GT. 0 .AND. NGROWD .LE. NGDAYS) THEN
            IF (DLAI(LDAY) .LT. PCD) THEN
              PTRANS(LDAY) = 0.0
            ELSE IF (DLAI(LDAY) .GT. PCE) THEN
              PTRANS(LDAY) = PET(LDAY)*(PCA+PCB*PCE**PCC)
            ELSE
              PTRANS(LDAY) = PET(LDAY)*(PCA+PCB*DLAI(LDAY)**PCC)
            ENDIF
            IF (PTRANS(LDAY) .GT. PET(LDAY)) PTRANS(LDAY) = PET(LDAY)
            PTRANS(LDAY) = (1.0-BARE)*PTRANS(LDAY)
            PEVAPO(LDAY) = PET(LDAY)-PTRANS(LDAY)
          ELSE
            PEVAPO(LDAY) = PET(LDAY)
          ENDIF
          TOTPET = TOTPET+PET(LDAY)
          TOTPTR = TOTPTR+PTRANS(LDAY)
          TOTPEV = TOTPEV+PEVAPO(LDAY)
        ENDDO
        IF (IPRINT .GT. 0) THEN
          IF (IPRINT .GT. 1) WRITE(LOUT,7220) (I,PET(I),                &
     &        PTRANS(I),PEVAPO(I),I=1,NDAYS)
          WRITE(LOUT,7230) TOTPET,TOTPTR,TOTPEV
        ENDIF
      ELSE IF (NFPET .EQ. 2) THEN
!
!     NFPET = 2  Daily PET and Hinds transpiration function to define  
!                PT and PE from PET (relation found by Ted Hinds).
!                In addition, the user can scale the data according
!                to biomass (BIOMAS).  The base biomass is 220 g/m**2
!
        IF (NSOW .GT. 90 .AND. NSOW .LT. 274) THEN
          IF (NSOW .LT. 152) THEN
            NSOW = 90
          ELSE 
            NSOW = 274
          ENDIF
          WRITE(LUS,7235) NSOW
          IF (LUSLUW) WRITE(LUS,7235) NSOW
        ENDIF
        IF (NHRVST .LT. 151 .OR. NHRVST .GT. 243) THEN
          NHRVST = 151
          WRITE(LUS,7236) NHRVST
          IF (LUSLUW) WRITE(LUS,7236) NHRVST
        ENDIF
        IF (BIOMAS .EQ. 0) BIOMAS = 220.0
        IF (IPRINT .GT. 1) THEN 
          WRITE(LOUT,7240) BIOMAS
          WRITE(LOUT,7210)
        ENDIF
        IF (NSOW .LE. 90) THEN
          NGDAYS = NHRVST-NSOW+1
          PREDAY = 90-NSOW
        ELSE 
          NGDAYS = ILEAP+365-NSOW+1+NHRVST
          PREDAY = ILEAP+365-NSOW+90+1
        ENDIF
        CRATIO = 0.300
        DO LDAY=1,NDAYS
          IF (NSOW .LE. LDAY) THEN
            NGROWD = LDAY-NSOW+1
          ELSE 
            NGROWD = ILEAP+365+LDAY-NSOW+1
          ENDIF
          IF (NGROWD .GT. 0 .AND. NGROWD .LE. NGDAYS) THEN
            IF (LDAY .LT. 90) THEN
              RATIO = CRATIO*NGROWD/PREDAY
            ELSE IF (LDAY .GE. 90 .AND. LDAY .LE. 151) THEN
              RATIO = CRATIO
            ELSE IF (LDAY .GT. 151 .AND. LDAY .LE. NHRVST) THEN
              RATIO = CRATIO*(NHRVST-LDAY+1)/(NHRVST-151+1)
            ELSE IF (LDAY .GT. 273) THEN
              RATIO = CRATIO*NGROWD/PREDAY
            ENDIF
            RATIO = RATIO*BIOMAS/220.0
            IF (RATIO .GT. 1.0) RATIO = 1.0
            PTRANS(LDAY) = (1.0-BARE)*RATIO*PET(LDAY)
            PEVAPO(LDAY) = PET(LDAY)-PTRANS(LDAY)
          ELSE
            PEVAPO(LDAY) = PET(LDAY)
          ENDIF
          TOTPET = TOTPET+PET(LDAY)
          TOTPTR = TOTPTR+PTRANS(LDAY)
          TOTPEV = TOTPEV+PEVAPO(LDAY)
        ENDDO
        IF (IPRINT .GT. 0) THEN
          IF (IPRINT .GT. 1) WRITE(LOUT,7220) (I,PET(I),                &
     &        PTRANS(I),PEVAPO(I),I=1,NDAYS)
          WRITE(LOUT,7230) TOTPET,TOTPTR,TOTPEV
        ENDIF
      ELSE 
        WRITE(LUS,7250) NFPET
        IF (LUSLUW) WRITE(LOUT,7250) NFPET
        GO TO 999 
      ENDIF
      RETURN
999   WRITE(LOUT,*) ' Program terminated prematurely in PETPART'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
7200  FORMAT(/,' NFPET = 1:',/,                                         &
     &'   PET is partitioned into PT and PE according to the',/,        &
     &'   relationship developed by Ritchie (1972)',/,                  &
     &'   The user-specified coefficients are:',/,                      &
     &'     a = ',F6.3,/,                                               &
     &'     b = ',F6.3,/,                                               &
     &'     c = ',F6.3,/,                                               &
     &'     d = ',F6.3,' (below this LAI, PT is zero)',/,               &
     &'     e = ',F6.3,' (above this LAI, PT=f(e))')
7210  FORMAT(/2(6X,'DAY',4X,'PET',3X,'PTRANS',2X,'PEVAPO'),/,           &
     &    2(6X,'---',2X,'------',2X,'------',2X,'------'))
7220  FORMAT(2(6X,I3,3F8.4))
7230  FORMAT(/,                                                         &
     &'   Totals:  PET    = ',F8.4,/,                                   &
     &'            PTRANS = ',F8.4,/,                                   &
     &'            PEVAPO = ',F8.4)
7235  FORMAT(' WARNING:  NSOW was out-of-bounds (91 to 273), so it',/,  &
     &' was reset to ',I3,/)
7236  FORMAT(' WARNING:  NHRVST was out-of-bounds (<151 or >243),',/,   &
     &       '           so it was reset to ',I3,/)
7240  FORMAT(' NFPET = 2:  PET is partitioned into PT and PE',/,        &
     &' according to the relationship developed by Ted Hinds',//,       &
     &' Biomass, g/m**2 (BIOMAS)    = ',F7.2,//)
7250  FORMAT(/,' ERROR:  Current NFPET options are 1 and 2.  Input',    &
     &' value was',I2,/)
!-----------------------------------------------------------------------
!     End of subroutine PETPART
!-----------------------------------------------------------------------
      END
      SUBROUTINE PRECIPBC(IPRINT,LUSLUW,LULI,LULO)
!-----------------------------------------------------------------------
!     Called by DATAINH and UNSATH
!
!     Read in precipitation data and write MET and precipitation data
!     to binary input file (LULO)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF(IPRINT .GT. 0)                                                 &
     &WRITE(LOUT,'(80("-"),/,A)')' Precipitation/irrigation parameters:'
!
!     If IRAIN = 0, precipitation data are provided; read NWATER
!
      IF (IRAIN .EQ. 0) THEN
        READ (LULI,'(I)') NWATER
        WRITE(LULO) NWATER
        IF (NWATER .EQ. 0) THEN
          IF (IPRINT .GT. 0)                                            &
     &      WRITE(LOUT,'(/,A)') '   NWATER = 0:  No rain or irrigation'
          IF (IETOPT .EQ. 0) GO TO 998
        ENDIF
      ELSE
!
!     NWATER is not needed for IRAIN > 0; written as dummy output
!
        WRITE(LULO) NWATER
      ENDIF
!
!     Water input:  IRDAY,IRTYP,NP,EFICEN,TIME,AMOUNT
!
      IF (IPRINT .GT. 1) THEN
        WRITE(LOUT,7409) IRAIN
        IF (IRAIN .EQ. 0) WRITE(LOUT,7410) NWATER
        WRITE(LOUT,7411)
      ENDIF
      IRDAY  = 0
      LIRDAY = IRDAY
      IRTYPE = 1
      EFICEN = 1.0
      DO I=1,NDAYS
        IF (IRAIN .EQ. 0) THEN
          IF ((IETOPT .EQ. 0 .AND. IRDAY .LT. I) .OR.                   &
     &      (IETOPT .EQ. 1 .AND. IRDAY .LT. RMDATA(1,I))) THEN
            READ(LULI,'(3I,D)',END=775) IRDAY,IRTYPE,NPP,EFICEN
            READ(LULI,'(2D)',END=775) (RTIME(J),AMOUNT(J),J=1,NPP)
          ENDIF
          IF (IETOPT .EQ. 0 .AND. IRDAY .LE. I) GO TO 777
          IF (IETOPT .EQ. 1 .AND. IRDAY .LE. RMDATA(1,I)) GO TO 777
775       NP = 0
          GO TO 769
777       NP = NPP
          IF (EFICEN .EQ. 0.0) EFICEN = 1.0
          DO J=1,NP-1
            TOTALR = TOTALR+AMOUNT(J)
            IF (RTIME(J+1) .LE. RTIME(J) .AND. IRTYPE .NE. 4) THEN
              WRITE(LUS,7430) IRDAY
              IF (LUSLUW) WRITE(LOUT,7430) IRDAY
              GO TO 999
            ENDIF
          ENDDO
          IF (IRDAY .LE. LIRDAY) THEN
            WRITE(LUS,7440) IRDAY
            IF (LUSLUW) WRITE(LOUT,7440) IRDAY
            GO TO 999
          ENDIF
          LIRDAY = IRDAY
        ELSE
          NP = 0
          RAINTOT = RMDATA(8,I)
          IF (RAINTOT .EQ. 0.0) GO TO 769
          NWATER = NWATER+1
          IRDAY = I
          IF (RAINTOT .LE. HPR) THEN
            NP = 2
            RTIME(1) = 0.0
            RTIME(2) = 1.0
            AMOUNT(1) = RAINTOT
            AMOUNT(2) = 0.0
          ELSE
            IF (24.0*HPR .lt. RAINTOT) THEN
              WRITE(*,7428) IRDAY,RAINTOT,HPR,HPR*24
7428  FORMAT(' Rain total exceeded 24*HPR on day',i4,/,                 &
     &' Rain total was',F8.4,/,                                         &
     &' HPR value and allowed amount = ',2F8.4,/,                       &
     &' For this day, HPR was set higher to allow full rain.')
              NP = 2
              RTIME(1) = 0.0
              RTIME(2) = 24.0
              AMOUNT(1) = RAINTOT
              AMOUNT(2) = 0.0
            ELSE
              FIRST = INT(RAINTOT/HPR)
              NP = 3
              RTIME(1) = 0.D0
              RTIME(2) = FIRST
              RTIME(3) = FIRST + 1.0
              AMOUNT(1) = FIRST*HPR
              AMOUNT(2) = RAINTOT-AMOUNT(1)
              AMOUNT(3) = 0.0
            ENDIF        
          ENDIF
          TOTALR = TOTALR+AMOUNT(1)+AMOUNT(2)
        ENDIF
        IF (IPRINT .GT. 1) THEN
          WRITE(LOUT,7450) IRDAY,RTIME(1),AMOUNT(1),IRTYPE,EFICEN,NP
          WRITE(LOUT,7460) (RTIME(J),AMOUNT(J),J=2,NP)
        ENDIF
769     IF (IETOPT .EQ. 0) THEN 
          IF (I .EQ. IRDAY) WRITE(LULO) IRDAY,IRTYPE,EFICEN,NP,         &
     &      (RTIME(J),AMOUNT(J),J=1,NP)
        ELSE
          METDAY = INT(RMDATA(1,I))
          WRITE(LULO) METDAY,(RMDATA(J,I),J=2,7),IRTYPE,EFICEN,NP,      &
     &      (RTIME(J),AMOUNT(J),J=1,NP)
        ENDIF
      ENDDO
      IF (IPRINT .GT. 0) THEN
        IF (IRAIN .EQ. 1) WRITE(LOUT,*)
        IF (IRAIN .EQ. 1) WRITE(LOUT,7410) NWATER
        WRITE(LOUT,7480) TOTALR
      ENDIF
998   RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in PRECIPBC'
      IF (LUSLUW) WRITE(LOUT,'(A)')                                     &
     &' Program terminated prematurely in PRECIPBC'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
7409  FORMAT(/,T4,'IRAIN =',I2,': precipitation data provided')
7410  FORMAT(/,T4,'NWATER (number of days of rain/irrigation) =',I3)
7411  FORMAT(/,T17,'Rainfall/Irrigation Details',//,                    &
     &T4,'Day',T11,'Time',T18,'Amount',T27,'Application',               &
     &T41,'Efficiency',T53,'Changes In',/,T11,'(hr)',T19,'(cm)',        &
     &T31,'Type',T54,'Rate/Head',/,T4,'---',T11,'----',T18,'------',    &
     &T29,'-----------',T41,'----------',T53,'----------')
7430  FORMAT(' ERROR:  Time of rain/irrigation not in ascending',       &
     &' order on day',I4)
7440  FORMAT(' ERROR:  Day of rain/irrigation not in ascending',        &
     &' order on day',I4)
7485  FORMAT(' ERROR:  IRTYPE cannot equal 4 while IETOPT',             &
     &' is greater than zero')
7450  FORMAT(T4,I3,T9,F6.3,T17,F7.4,T32,I1,T43,F5.3,T57,I2)
7460  FORMAT(T9,F6.3,T17,F7.4)
7480  FORMAT(/,T4,'Total Water Applied (cm) = ',1PE13.6)
!-----------------------------------------------------------------------
!     End of subroutine PRECIPBC
!-----------------------------------------------------------------------
      END
      SUBROUTINE HEATPAR(IPRINT,LUSLUW,MATN)
!-----------------------------------------------------------------------
!     called by DATAINH
!
!     reads the soil parameters for heat flow
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (IPRINT .GT. 0) WRITE(LOUT,4000)
      DO J=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(6D)',ERR=999) (TCON(I,J),I=1,5),CHS(J)
        IF (IPRINT .GT. 0)                                              &
     &    WRITE(LOUT,4350) DUMMY,J,(TCON(I,J),I=1,5),CHS(J)
!
!     The thermal conductivity and enhancement factor parameters are 
!     reduced. The conductivity units are changed from W/m/K to J/cm/h/K
!
        THETJ = THET(J)
        TCON(4,J) = 3.6E3*(TCON(1,J)-TCON(4,J))/1.E2
        TCON(1,J) = 3.6E3*TCON(1,J)/1.E2
        TCON(2,J) = 3.6E3*(TCON(2,J)/THETJ)/1.E2
        TCON(3,J) = TCON(3,J)/THETJ
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(5D)',ERR=999) (EF(I,J),I=1,5)
        IF (IPRINT .GT. 0) WRITE(LOUT,4350) DUMMY,J,(EF(I,J),I=1,5)
        THETJ = THET(J)
        EF(2,J) = EF(2,J)/THETJ
        EF(4,J) = EF(1,J)-EF(4,J)
        EF(3,J) = EF(3,J)/THETJ
      ENDDO
      DO I=1,NPT
        MATI = MAT(I)
        CHSOIL(I) = CHS(MATI)*(1-THET(MATI))
      ENDDO
      RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in HEATPAR'
      IF (LUSLUW) WRITE(LOUT,'(A)')                                     &
     &' Program terminated prematurely in HEATPAR'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
4000  FORMAT(80("-"),/,' Soil thermal properties:',//,
     &'   IHEAT = 1:  Thermal conductivity and heat capacity',          &
     &' parameters',/,T64,'Heat',/,                                     &
     &' Material',T15,'A',T25,'B',T35,'C',T45,'D',T55,'E',T62,          &
     &'Capacity',/,T62,'(J/cm3/K)')                                     &
4350  FORMAT(1X,A60,/,' No.',I5,1X,T11,1P7E10.3)
!-----------------------------------------------------------------------
!     End of subroutine HEATPAR
!-----------------------------------------------------------------------
      END
      SUBROUTINE HYSPAR(IPRINT,LUSLUW,MATN)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     Reads soil parameters for hysteresis
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (IPRINT .GT. 0) WRITE(LOUT,4351)
      DO M=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(I,3D)',ERR=999)IPATHA(M),SARWA(M),ALFACT(M),HYSHPH(M)
        IF (IPATHA(M) .EQ. 2 .OR. IPATHA(M) .GT. 7) THEN
          WRITE(LUS,'(A)') ' IPATHA can only be 1, 3, 4, 5, 6, or 7'
          GOTO 999
        ENDIF
        IF (ALFACT(M) .LT. 1.0) THEN
          ALFACT(M) = 1.0
          WRITE(LUS,'(A)') ' ALFACT reset to 1.0 for material ',M
        ENDIF
        IF (HYSHPH(M) .LT. 0.0 .OR. HYSHPH(M) .GT. HYSMXH) THEN
          HYSHPH(M) = HYSMXH
          WRITE(LUS,'(A)') ' HYSHPH reset to HYSMXH for material ',M
        ENDIF
        IF (IPRINT .GT. 0) WRITE(LOUT,4352) DUMMY,M,IPATHA(M),          &
     &    SARWA(M),ALFACT(M),HYSHPH(M)
!
!     Transform SARWA from saturation to effective saturation
!
        SARWA(M) = SARWA(M)/(1.0-THTR(M)/THET(M))
      ENDDO
      RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in HYSPAR'
      IF (LUSLUW)                                                       &
     &WRITE(LOUT,'(A)') ' Program terminated prematurely in HYSPAR'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
4351  FORMAT(80("-"),/,' Hysteresis parameters:',//,'   IHYS > 0',/,    &
     &T16,'MAX NO.',/,'   Material',T16,'PATHS',T30,'SARWA',            &
     &T45,'ALFACT',T60,'HYSHPH')
4352  FORMAT(3X,A77,/,'   No. ',I2,T17,I2,T24,1P3E15.7)
!-----------------------------------------------------------------------
!     End of subroutine HYSPAR
!-----------------------------------------------------------------------
      END
      SUBROUTINE SHPPAR(IPRINT,LUSLUW,MATN)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     Read in soil hydraulic property parameters
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      FUNWH(HTMP,VGA,VGN) = 1./(1.+(VGA*HTMP)**VGN)
      INC = M4
      IF (IPRINT .GT. 0)                                                &
     &  WRITE(LOUT,'(/,A,/)') ' Soil hydraulic properties:'
      IF (KOPT .NE. 1) GOTO 400
!-----------------------------------------------------------------------
!     KOPT = 1: Polynomial functions
!-----------------------------------------------------------------------
      READ(LUI,'(2I)',ERR=999) MAXPOL,MAXCOE
      IF (IPRINT .GT. 0) WRITE(LOUT,4000) MAXPOL,MAXCOE
      IF (MAXPOL .GT. M3) THEN
        WRITE(LUS,'(A)') ' MAXPOL > M3 dimension in unsath.inc file'
        GOTO 999
      ENDIF
      INC = MAXPOL*MAXCOE
      IF (INC .GT. M4) THEN
        WRITE(LUS,'(A)') ' MAXPOL*MAXCOE > M4 dim. in unsath.inc file'
        GOTO 999
      ENDIF
      DO I=1 ,MATN
      READ(LUI,'(A)',ERR=999) DUMMY
      READ (LUI,'(I,2D)',ERR=999)NSUBTH(I),AIRINT(I),THET(I)
      IF (IPRINT .GT. 0) WRITE(LOUT,4020) I,DUMMY,NSUBTH(I),AIRINT(I),  &
     &  THET(I)
      DO J=1,NSUBTH(I)
        ID = (I-1)*MAXPOL+J
        IP = (I-1)*INC+(J-1)*MAXCOE
        READ(LUI,'(2I,2D)',ERR=999) II,NDEGTH(ID),XX,XDIVTH(ID)
        IF (J .GT. 1) THEN
          IF (XDIVTH(ID) .LT. XDIVTH(ID-1)) THEN
            WRITE(LUS,4040) I,J,XDIVTH(ID),XDIVTH(ID-1)
            IF (LUSLUW) WRITE(LOUT,4040) I,J,XDIVTH(ID),XDIVTH(ID-1)
            GOTO 999
          ENDIF
        ENDIF
        READ(LUI,'(5D)',ERR=999) (CREGTH(IP+L),L=1,NDEGTH(ID))
        IF (IPRINT .GT. 0) WRITE(LOUT,4060) ID,NDEGTH(ID),XDIVTH(ID),   &
     &    (CREGTH(IP+L),L=1,NDEGTH(ID))
      ENDDO  !  Loop_J
4000  FORMAT('   KOPT = 1:  Polynomial hydraulic functions',            &
     &T15,'MAXPOL = ',I1,', MAXCOE = ',I1)
4020  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T4,'NSUBTH = ',I1,', AIRINT = ',F9.4,', THET = ',F6.4,//,         &
     &T4,'ID  NDEGTH     XDIVTH        CREGTH',/,                       &
     &T4,'--  ------  ------------  ------------')
4060  FORMAT(0P,T3,I2,T9,I2,T15,1P,E12.4,(T29,E13.6))
4040  FORMAT(' ERROR:  Material No.',I3,', J =',I4,', XDIVTH(J) = ',    &
     &E15.7,' < XDIVTH(J-1) = ',E15.7)
      READ(LUI,'(A)',ERR=999) DUMMY
      READ(LUI,'(I,2D)',ERR=999) NSUBKH(I),AIRINK(I),SK(I)
      IF (IPRINT .GT. 0)                                                &
     &    WRITE(LOUT,4070) DUMMY,NSUBKH(I),AIRINK(I),SK(I)
      DO J=1,NSUBKH(I)
        ID = (I-1)*MAXPOL+J
        IP = (I-1)*INC+(J-1)*MAXCOE
        READ(LUI,'(2I,2D)',ERR=999) II,NDEGKH(ID),XX,XDIVKH(ID)
        IF (J .GT. 1) THEN
          IF (XDIVKH(ID) .LT. XDIVKH(ID-1)) THEN
            WRITE(LUS,4080) I,J,XDIVKH(ID),XDIVKH(ID-1)
            IF (LUSLUW) WRITE(LOUT,4080) I,J,XDIVKH(ID),XDIVKH(ID-1)
            GOTO 999
          ENDIF
        ENDIF
        READ(LUI,'(5D)',ERR=999) (CREGKH(IP+L),L=1,NDEGKH(ID))
        IF (IPRINT .GT. 0) WRITE(LOUT,4060) ID,NDEGKH(ID),XDIVKH(ID),   &
     &    (CREGKH(IP+L),L=1,NDEGKH(ID))
      ENDDO
      ENDDO
4080  FORMAT(' ERROR:  Material No. ',I2,', J = ',I3,', XDIVKH(J) = ',  &
     &E15.7,' < XDIVKH(J-1) = ',E15.7)
4070  FORMAT('   K = f(H), ',A60,1P,//,                                 &
     &T4,'NSUBKH = ',I1,', AIRINK = ',F9.4,', SK = ',1P,E12.4,/,        &
     &T4,'ID  NDEGKH     XDIVKH        CREGKH',/,                       &
     &T4,'--  ------  ------------  ------------')
      GOTO 440
400   IF (KOPT .NE. 2) GOTO 450
!-----------------------------------------------------------------------
!     KOPT = 2: Haverkamp hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4100)
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(6D)',ERR=999)                                        &
     &    THET(I),THTR(I),AIRINT(I),ALPHA,BETA,RETOPT
        IF (IPRINT .GT. 0) WRITE(LOUT,4110) I,DUMMY,THET(I),THTR(I),    &
     &    AIRINT(I),ALPHA,BETA,RETOPT
        IP = (I-1)*INC
        CREGTH(IP+1) = ALPHA
        CREGTH(IP+2) = BETA
        CREGTH(IP+3) = RETOPT
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(4D)',ERR=999) AIRINK(I),SK(I),A,B
        IF (IPRINT .GT. 0) WRITE(LOUT,4120) DUMMY,AIRINK(I),SK(I),A,B
        CREGKH(IP+1) = A
        CREGKH(IP+2) = B
      ENDDO
4100  FORMAT('   KOPT = 2: Haverkamp hydraulic functions')
4110  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T14,'  THET =',E12.5,T36,'  THTR =',E12.5,T58,'AIRINT =',E12.5,/, &
     &T14,' ALPHA =',E12.5,T36,'  BETA =',E12.5,T58,'RETOPT =',E12.5)
4120  FORMAT('   K = f(H), ',A60,1P,/,                                  &
     &T14,'AIRINK =',E12.5,T36,'    SK =',E12.5,T58,'     A =',E12.5,/, &
     &T14,'     B =',E12.5)
      GOTO 440
450   IF (KOPT .NE. 3) GOTO 470
!-----------------------------------------------------------------------
!     KOPT = 3: Brooks-Corey hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4200)
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(4D)',ERR=999) THET(I),THTR(I),AIRINT(I),B
        IF (IPRINT .GT. 0)                                              &
     &      WRITE(LOUT,4210) I,DUMMY,THET(I),THTR(I),AIRINT(I),B
        IP = (I-1)*INC
        CREGTH(IP+1) = B
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(5D)',ERR=999) RKMOD,SK(I),AIRINK(I),B,EPIT
        IF (IPRINT .GT. 0)                                              &
     &      WRITE(LOUT,4220) DUMMY,AIRINK(I),B,SK(I),EPIT,RKMOD
        CREGKH(IP+1) = RKMOD
        CREGKH(IP+2) = B
        CREGKH(IP+3) = EPIT
      ENDDO
4200  FORMAT('   KOPT = 3: Brooks-Corey hydraulic functions')
4210  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T14,'  THET =',E12.5,T36,'  THTR =',E12.5,/,                      &
     &T14,'AIRINT =',E12.5,T36,'     B =',E12.5)
4220  FORMAT('   K = f(H), ',A60,1P,/,                                  &
     &T14,'AIRINK =',E12.5,T36,'     B =',E12.5,T58,'    SK =',E12.5,/, &
     &T14,'  EPIT =',E12.5,T36,' RKMOD =',E12.5)
      GOTO 440
470   IF (KOPT .NE. 4) GOTO 487
!-----------------------------------------------------------------------
!     KOPT = 4: van Genuchten hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4222)
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(4D)',ERR=999) THET(I),THTR(I),VGA,VGN
        AIRINT(I) = 0.0
        VGM = 1.0-1.0/VGN
        IF (IPRINT .GT. 0) WRITE(LOUT,4224) I,DUMMY,                    &
     &      THET(I),THTR(I),VGA,VGN,VGM
        IP = (I-1)*INC
        CREGTH(IP+1) = VGA
        CREGTH(IP+2) = VGN
        CREGTH(IP+3) = VGM
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(5D)',ERR=999) RKMOD,SK(I),VGA,VGN,EPIT
        AIRINK(I) = 0.0
        VGM = 1.0-1.0/VGN
        IF (IPRINT .GT. 0) WRITE(LOUT,4226) DUMMY,RKMOD,SK(I),VGA,VGN,  &
     &    VGM,EPIT
        CREGKH(IP+1) = VGA
        CREGKH(IP+2) = VGN
        CREGKH(IP+3) = VGM
        CREGKH(IP+4) = RKMOD
        CREGKH(IP+5) = EPIT
      ENDDO
4222  FORMAT('   KOPT = 4: van Genuchten hydraulic functions')
4224  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T14,'  THET =',G12.5,T36,'  THTR =',G12.5,T58,' ALPHA =',G12.5,/, &
     &T14,'     N =',G12.5,T36,'     M =',G12.5)
4226  FORMAT('   K =f(H), ',A60,/,1P,                                   &
     &T14,' RKMOD =',G12.5,T36,'    SK =',G12.5,T58,'     A =',G12.5,/, &
     &T14,'     N =',G12.5,T36,'     M =',G12.5,T58,'  EPIT =',G12.5)
      GOTO 440
487   IF (KOPT .NE. 5) GOTO 488
!-----------------------------------------------------------------------
!     KOPT = 5: Fayer-Simmons modified Brooks-Corey hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4212)
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(5D)',ERR=999) THET(I),THTA(I),AIRINT(I),B,HM(I)
        IF (IPRINT .GT. 0)                                              &
     &      WRITE(LOUT,4214) I,DUMMY,AIRINT(I),B,THET(I),THTA(I),HM(I)
        IP = (I-1)*INC
        CREGTH(IP+1) = B
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(2D)',ERR=999) SK(I),EPIT
        IF (IPRINT .GT. 0)                                              &
     &      WRITE(LOUT,4216) DUMMY,SK(I),EPIT
        SA  = THTA(I)/THET(I)
        HE  = AIRINT(I)
        HMM = HM(I)
        A   = HE**(1/B)
        P   = -1/B
        PM1 = P-1
        BC1 = SA/LOG(HMM)
        BC2 = A*P*BC1
        BC3 = A*(P-P*SA+BC1)-BC2/PM1
        GAMMAM = (BC3+BC2*LOG(HMM))*HMM**PM1/PM1+BC1/HMM
        GAMMAX = (BC3+BC2*LOG( HE))* HE**PM1/PM1+BC1/HE-GAMMAM
        CREGKH(IP+1) = EPIT
        CREGKH(IP+2) = BC1
        CREGKH(IP+3) = BC2
        CREGKH(IP+4) = BC3
        CREGKH(IP+5) = PM1
        CREGKH(IP+6) = GAMMAM
        CREGKH(IP+7) = GAMMAX
      ENDDO
4212  FORMAT('   KOPT =5:',                                             &
     &'  Fayer-Simmons modified Brooks-Corey hydraulic functions')
4214  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T14,'  THET =',G12.5,T36,'  THTA =',G12.5,/,T58,'AIRINT =',G12.5, &
     &T14,'     B =',G12.5,T36,'    HM =',G12.5)
4216  FORMAT('   K =f(H), ',A60,/,1P,
     &T14,'    SK =',G12.5,T36,'  EPIT =',G12.5)
      GOTO 440
488   IF (KOPT .NE. 6) GOTO 491
!-----------------------------------------------------------------------
!     KOPT = 6: Fayer-Simmons modified van Genuchten hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4223)
      F  = 0.01
      WO = 1.E-10
      XO = 0.9
      IF (IPRINT .GT. 0) WRITE(LOUT,4229) F,WO,XO
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(5D)',ERR=999) THET(I),THTA(I),VGA,VGN,HM(I)
        VGM = 1.0-1.0/VGN
        IF (IPRINT .GT. 0)                                              &
     &    WRITE(LOUT,4225) I,DUMMY,THET(I),THTA(I),VGA,VGN,VGM,HM(I)
        IP = (I-1)*INC
        HMM = HM(I)
        AMN = VGA*VGM*VGN
        RLOGHM = LOG(HMM)
        THTAHM = THTA(I)/RLOGHM
        CREGTH(IP+1)  = VGA
        CREGTH(IP+2)  = VGN
        CREGTH(IP+3)  = VGM
        CREGTH(IP+6)  = AMN
        CREGTH(IP+7)  = RLOGHM
        CREGTH(IP+8)  = THTAHM
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(2D)',ERR=999) SK(I),EPIT
        IF (IPRINT .GT. 0) WRITE(LOUT,4227) DUMMY,SK(I),EPIT
        SA    = THTA(I)/THET(I)
        WHM   = FUNWH(HMM,VGA,VGN)
        RLAMB = 1.0/LOG(HMM)
        GAMMA = RLAMB*LOG(VGA)
        GAMSA = GAMMA*SA
!
!     Calculate GAMMAM
!
        GAMMAM = VGA*(VGN-1)*(VGA*HMM)**(-VGN)/VGN
!
!     Calculate HO and  HC
!
        HO = F/VGA
        XC = (1.-GAMSA+RLAMB*SA/VGN)/SA
        HC = HMM**(1.-XC)/VGA
        IF (HC .LT. 1.E-20) HC = 1.E-20
        WHO = FUNWH(HO,VGA,VGN)
        WHC = FUNWH(HC,VGA,VGN)
        AHONM1 = (VGA*HO)**(VGN-1.)
        AHCNM1 = (VGA*HC)**(VGN-1.)
        RLAHON = LOG(VGA*HO)-1./(VGN-1.)
        RLAHCN = LOG(VGA*HC)-1./(VGN-1.)
        GAMMAO = VGA*((1.-GAMSA-SA+RLAMB*SA/VGN)*(AHONM1-AHCNM1)        &
     &           +RLAMB*SA*(AHONM1*RLAHON-AHCNM1*RLAHCN))
        P   = VGM+1.
        CK1 = -VGM
        CK2 = -CK1*(VGM-1.)/2.
        CK3 = -CK2*(VGM-2.)/3.
        IF (WHM .GT. WO) WO = WHM
        CALL FHX(1.-XO,HXO)
        CALL FGSX(XO,GSXO)
        RI3BWO = ((WO-WHM)+WHM*LOG(WHM)-WO*LOG(WO))*VGM
        CALL INTGRL(HO,WHO,GAMMAS)
        GAMMAX = GAMMAS + GAMMAO
        CREGKH(IP+1)  = EPIT
        CREGKH(IP+2)  = WO
        CREGKH(IP+3)  = GAMMAM
        CREGKH(IP+4)  = GAMMA
        CREGKH(IP+5)  = RLAMB
        CREGKH(IP+6)  = SA
        CREGKH(IP+7)  = WHM
        CREGKH(IP+8)  = CK1
        CREGKH(IP+9)  = CK2
        CREGKH(IP+10) = CK3
        CREGKH(IP+11) = HXO
        CREGKH(IP+12) = GSXO
        CREGKH(IP+13) = RI3BWO
        CREGKH(IP+14) = GAMMAX
      ENDDO
4223  FORMAT('   KOPT = 6: modified van Genuchten hydraulic functions')
4225  FORMAT(/,'   Material No.',I5,/,'   THETA =f(H), ',A60,1P,/,      &
     &T14,'  THET =',E12.5,T36,'  THTA =',E12.5,T58,' ALPHA =',E12.5,/, &
     &T14,'     N =',E12.5,T36,'     M =',E12.5,T58,'    HM =',E12.5)
4227  FORMAT('   K =f(H), ',A60,/,1P,                                   &
     &T14,'    SK =',E12.5,T36,'  EPIT =',E12.5)
4229  FORMAT(/,'   Parameters F, WO, and XO are used for KOPT = 6',/,   &
     &T4,'calculations of conductivity.  For all materials,',/,         &
     &T4,'the values are constant.  They are:',/,                       &
     &T14,'     F =',E12.5,T36,'    WO =',E12.5,T58,'    XO =',E12.5)
      GOTO 440
491   IF (KOPT .NE. 7) GOTO 490
!-----------------------------------------------------------------------
!     KOPT = 7: Rossi-Nimmo "sum" hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4232)
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(5D)',ERR=999) THET(I),PSID,PSIO,RLAM,PSII
        HOHIL = (PSIO/PSII)**RLAM
        RNA = -(0.5*RLAM*HOHIL-1+HOHIL-(PSIO/PSID)**RLAM)               &
     &        /(0.5+LOG(PSID/PSII))
        RNC = 0.5*(RLAM*HOHIL+RNA)*(PSIO/PSII)**2
        IF (IPRINT .GT. 0)                                              &
     &    WRITE(LOUT,4234) I,DUMMY,THET(I),PSID,PSIO,RLAM,RNA,RNC,PSII
        IP = (I-1)*INC
        CREGTH(IP+1) = RNA
        CREGTH(IP+2) = RNC
        CREGTH(IP+3) = PSID
        CREGTH(IP+4) = PSIO
        CREGTH(IP+5) = RLAM
        CREGTH(IP+6) = PSII
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(2D)',ERR=999) SK(I),EPIT
        IF (IPRINT .GT. 0) WRITE(LOUT,4235) DUMMY,SK(I),EPIT
        G2CON = (RNA+(RLAM/(1+RLAM))*(PSIO/PSID)**RLAM)/PSID
        GAMMA2 = (RNA+(RLAM/(1+RLAM))*HOHIL)/PSII-G2CON
        GAMMAX = GAMMA2+2*RNC*PSII/PSIO**2
        CREGKH(IP+1) = EPIT
        CREGKH(IP+2) = G2CON
        CREGKH(IP+3) = GAMMA2
        CREGKH(IP+4) = GAMMAX
      ENDDO
4232  FORMAT('   KOPT = 7: Rossi-Nimmo "sum" hydraulic functions')
4234  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T14,'  THET =',E12.5,T36,'  PSID =',E12.5,T58,'  PSIO =',E12.5,/, &
     &T14,'  RLAM =',E12.5,T36,'   RNA =',E12.5,T58,'   RNC =',E12.5,/, &
     &T14,'  PSII =',E12.5)
4235  FORMAT('   K =f(H), ',A60,/,1P,                                   &
     &T14,'    SK =',E12.5,T36,'  EPIT =',E12.5)
      GOTO 440
490   IF (KOPT .NE. 8) THEN
        WRITE(LUS,4230) KOPT
        IF (LUSLUW) WRITE(LOUT,4230) KOPT
4230  FORMAT(' ERROR:  Current KOPT options are 1 to 8.',/,             &
     &       '         The input value was ',I3)
        GOTO 999
      ENDIF
!-----------------------------------------------------------------------
!     KOPT = 8: Rossi-Nimmo "junction" hydraulic functions
!-----------------------------------------------------------------------
      IF (IPRINT .GT. 0) WRITE(LOUT,4236)
      DO I=1,MATN
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(4D)',ERR=999) THET(I),PSID,PSIO,RLAM
        PSII = PSIO*(0.5*RLAM+1)**(1/RLAM)
        PSIJ = PSID*EXP(-1/RLAM)
        RNA  = (PSIO/PSIJ)**RLAM/LOG(PSID/PSIJ)
        RNC  = (1-(PSIO/PSII)**RLAM)*(PSIO/PSII)**2
        IF (IPRINT .GT. 0) WRITE(LOUT,4237)                             &
     &     I,DUMMY,THET(I),PSID,PSIO,RLAM,RNA,RNC,PSII,PSIJ
        IP = (I-1)*INC
        CREGTH(IP+1) = RNA
        CREGTH(IP+2) = RNC
        CREGTH(IP+3) = PSID
        CREGTH(IP+4) = PSIO
        CREGTH(IP+5) = RLAM
        CREGTH(IP+6) = PSII
        CREGTH(IP+7) = PSIJ
        READ(LUI,'(A)',ERR=999) DUMMY
        READ(LUI,'(2D)',ERR=999) SK(I),EPIT
        IF (IPRINT .GT. 0) WRITE(LOUT,4238) DUMMY,SK(I),EPIT
        BETA   = (RLAM/(1+RLAM))*PSIO**RLAM
        GAMMA3 = RNA*((1.0/PSIJ)-(1.0/PSID))
        GAMMA2 = BETA*(PSII**(-RLAM-1)-PSIJ**(-RLAM-1))
        GAMMAX = GAMMA3+GAMMA2+2*RNC*PSII/PSIO**2
        CREGKH(IP+1) = EPIT
        CREGKH(IP+2) = BETA
        CREGKH(IP+3) = GAMMA3
        CREGKH(IP+4) = GAMMA2
        CREGKH(IP+5) = GAMMAX
      ENDDO
4236  FORMAT('   KOPT = 8: Rossi-Nimmo "junction" hydraulic functions')
4237  FORMAT(/,'   Material No.',I5,/,'   THETA = f(H), ',A60,1P,/,     &
     &T14,'  THET =',E12.5,T36,'  PSID =',E12.5,T58,'  PSIO =',E12.5,/, &
     &T14,'  RLAM =',E12.5,T36,'   RNA =',E12.5,T58,'   RNC =',E12.5,/, &
     &T14,'  PSII =',E12.5,T36,'  PSIJ =',E12.5)
4238  FORMAT('   K = f(H), ',A60,/,1P,                                  &
     &T14,'    SK =',E12.5,T36,'  EPIT =',E12.5)
440   RETURN
999   WRITE(LUS,'(A)') ' Program terminated prematurely in SHPPAR'
      IF (LUSLUW) WRITE(LOUT,'(A)')                                     &
     &' Program terminated prematurely in SHPPAR'
      STOP
!-----------------------------------------------------------------------
!     End of subroutine SHPPAR
!-----------------------------------------------------------------------
      END
      SUBROUTINE FILBLDR(IFLIST,NFLIST,IPRINT,ISTEAD)
!-----------------------------------------------------------------------
!     Called by DATAINH
!
!     reads in list of files containing P(ET)/M(ET) and PRECIP data
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      CHARACTER*80 DIRLOC,FILPRE,FILEXT,CIYS,FNAMEL
      IF (IFLIST .EQ. 1) THEN
        NY = NYEARS
        IF (ISTEAD .EQ. 1) NY = 1
      ELSE IF (IFLIST .EQ. 2 .OR. IFLIST .EQ. 3) THEN
        NY = NFLIST
      ENDIF
      IF (IFLIST .EQ. 1) THEN
!
!     Read the file descriptors for the PET/MET files
!
        READ(LUI,'(A)',ERR=999) DIRLOC,FILPRE,FILEXT
        I = 0
100     IF (I .EQ. NY) GO TO 200
        IYEAR = IYS+I
        I = I+1
        WRITE(UNIT=CIYS,FMT='(I)') IYEAR
        CALL CLIP(CIYS)
        ICHAR = INDEX(DIRLOC," ")
        PMFN(I)(1:ICHAR) = DIRLOC
        JCHAR = INDEX(FILPRE," ")
        PMFN(I)(ICHAR:ICHAR+JCHAR-1) = FILPRE
        JCHAR = INDEX(CIYS," ")
        KCHAR = INDEX(PMFN(I)," ")
        PMFN(I)(KCHAR:KCHAR+JCHAR-1) = CIYS
        KCHAR = INDEX(PMFN(I)," ")
        PMFN(I)(KCHAR:KCHAR) = "."
        JCHAR = INDEX(FILEXT," ")
        PMFN(I)(KCHAR+1:KCHAR+JCHAR-1) = FILEXT
        GO TO 100
!
!     Read the file descriptors for the PRECIP files
!
200     IF (IRAIN .EQ. 1) GOTO 400
        READ(LUI,'(A)',ERR=999) DIRLOC,FILPRE,FILEXT
        I = 0
300     IYEAR = IYS+I
        I = I+1
        WRITE(UNIT=CIYS,FMT='(I)') IYEAR
        CALL CLIP(CIYS)
        ICHAR = INDEX(DIRLOC," ")
        PRFN(I)(1:ICHAR) = DIRLOC
        JCHAR = INDEX(FILPRE," ")
        PRFN(I)(ICHAR:ICHAR+JCHAR-1) = FILPRE
        JCHAR = INDEX(CIYS," ")
        KCHAR = INDEX(PRFN(I)," ")
        PRFN(I)(KCHAR:KCHAR+JCHAR-1) = CIYS
        KCHAR = INDEX(PRFN(I)," ")
        PRFN(I)(KCHAR:KCHAR) = "."
        JCHAR = INDEX(FILEXT," ")
        PRFN(I)(KCHAR+1:KCHAR+JCHAR-1) = FILEXT
        IF (I .LT. NY) GO TO 300
      ELSE IF (IFLIST .EQ. 2) THEN
!
!     Read the file names from the input file
!
        READ (LUI,'(A)',ERR=999) (PMFN(I),I=1,NY)
        IF (IRAIN .EQ. 0) READ (LUI,'(A)',ERR=999) (PRFN(I),I=1,NY)
      ELSE IF (IFLIST .EQ. 3) THEN
!
!     Read the file names from external files (named in the input file)
!
        READ(LUI,'(A)',ERR=999) FNAMEL
        OPEN(UNIT=LUX,FILE=FNAMEL,STATUS='OLD',IOSTAT=IOS)
        IF (IOS .NE. 0) THEN
          WRITE(LUS,*) ' Unable to open file ',FNAMEL
          GO TO 999
        ENDIF
        READ (LUX,'(A)',ERR=999) (PMFN(I),I=1,NY)
        IF (IRAIN .EQ. 0) READ (LUX,'(A)',ERR=999) (PRFN(I),I=1,NY)
        CLOSE (UNIT=LUX)
      ENDIF
400   IBADFIL = 0
      DO I=1,NY
        OPEN(UNIT=LUX,FILE=PMFN(I),STATUS='OLD',IOSTAT=IOS)
        IF (IOS .NE. 0) THEN
          WRITE(LUS,6020) PMFN(I)
          IF (LUSLUW) WRITE(LUW,6020) PMFN(I)
          IBADFIL = IBADFIL+1
        ENDIF
        CLOSE (UNIT=LUX)
        IF (IRAIN .EQ. 0) THEN
          OPEN(UNIT=LUX,FILE=PRFN(I),STATUS='OLD',IOSTAT=IOS)
          IF (IOS .NE. 0) THEN
            WRITE(LUS,6025) PRFN(I)
            IF (LUSLUW) WRITE(LUW,6025) PRFN(I)
            IBADFIL = IBADFIL+1
          ENDIF
          CLOSE (UNIT=LUX)
        ENDIF
      ENDDO
      IF (IBADFIL .GT. 0) THEN
        WRITE(LUS,6030) IBADFIL
        IF (LUSLUW) WRITE(LUW,6030) IBADFIL
        GO TO 999
      ENDIF
      IF (IPRINT .GT. 1) THEN
        WRITE(LOUT,6000) (I,PMFN(I),I=1,NY)
        IF (IRAIN .EQ. 0) WRITE(LOUT,6010) (I,PRFN(I),I=1,NY)
      ENDIF
      RETURN
999   WRITE(LUS,'(A)') ' Subroutine FILBLDR terminated prematurely'
      IF (LUSLUW)                                                       &
     &  WRITE(LOUT,'(A)') ' Subroutine FILBLDR terminated prematurely'
      STOP
!-----------------------------------------------------------------------
!     Format Statements
!-----------------------------------------------------------------------
6000  FORMAT(80("-"),//' Multiyear Option',/,                           &
     &' Sim Year    PET input file name',/,(I9,4X,A65))
6010  FORMAT(' Multiyear Option',/,                                     &
     &' Sim Year    Precipitation input file name',/,(I9,4X,A65))
6020  FORMAT(' Could not read PET/MET file ',A50)
6025  FORMAT(' Could not read precip. file ',A50)
6030  FORMAT(' Total number of files that were unreadable =',I5)
!-----------------------------------------------------------------------
!     End of FILBLDR
!-----------------------------------------------------------------------
      END
      SUBROUTINE CLIP(ACHR)
!-----------------------------------------------------------------------
!     Removes unwanted characters in file names
!-----------------------------------------------------------------------
      CHARACTER*80 ACHR,TEMP
      ISTART = 1
      INDICE = 0
100   IF (ICHAR(ACHR(ISTART:ISTART)) .EQ. 32) THEN
        ISTART = ISTART+1
      ELSE
        INDICE = INDICE+1
        TEMP(INDICE:INDICE) = ACHR(ISTART:ISTART)
        ISTART = ISTART+1
      ENDIF
      IF (ISTART .LE. 80) GO TO 100
      ACHR(1:INDICE) = TEMP(1:INDICE)
      ACHR(INDICE+1:80) = " "
      RETURN
!-----------------------------------------------------------------------
!     End of CLIP
!-----------------------------------------------------------------------
      END
