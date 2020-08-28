      PROGRAM UNSATH
!-----------------------------------------------------------------------
!     Calls AIRTMP, BLR, DELCHK, FLUX, HEATFLOW, HDRYCALC, INTERK, 
!     KVCALC, POLYKH, RESET, RETENT, RLD, THERME, TIMEX, TRIDAG,
!     ROOTS, TRSINK, ZEROA, ZEROI, ZEROR, RELHUM
!
!     Solves water and heat balance equations 
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'      
      INTEGER SUBSTP,SUBUBC,WINDEX,SUBAST,SUBSDH
      CHARACTER*80 AYEAR
      LOGICAL FINISH
      DIMENSION A1(M1),A2(M1),A3(M1),ZY(M1),RATE(24),HSOL(M1),SUBQL(M1),&
     &    SUBQVH(M1),SUBQVT(M1),SINK(M1),SUBSNK(M1)
      FUNVD(TD) = EXP(46.440973-(6790.4985/TD)-6.02808*LOG(TD))
      DATA UV/3.01/,HLB/-20000.0/
      INCLUDE 'init.inc'
!
      WRITE(LUS,9000)
9000  FORMAT(/' Results are based on the use of unverified software.'/, &
     &    ' No assurance is expressed or implied as to the accuracy,'/, &
     &    ' completeness or usefulness of this information.'/)
      FINISH = .FALSE.
      CALL TIMEX
!
!     BIFILE = Binary input filename
!
      CALL WELCOME('UNSAT-H',UV)
      WRITE(LUS,8) 
8     FORMAT(/,' Enter input filename without the "*.BIN" extension',/, &
     &' (a "0" terminates the program) ===> ',$)
      READ(LUR,'(A80)',END=999) TMFILE
      IF(INDEX(TMFILE(1:1),'0') .GT. 0) GO TO 999
      RFILE = TMFILE
      NCHR = INDEX(TMFILE,' ')
      TMFILE(NCHR:NCHR+3) = '.bin'
      OPEN(UNIT=LUI,FILE=TMFILE,STATUS='OLD',FORM='UNFORMATTED')
      INQUIRE(UNIT=LUI,NAME=BIFILE)
      WRITE(LUS,10) SDATE,STIME,BIFILE
10    FORMAT(/,
     &' Date: ',A11,/,
     &' Time: ',A11,/,
     &' Binary Input Filename:',/,1X,A79,/)
      WRITE(*,*) ' Input for screen output (0=none; 1=yearly; 2=daily)'
      READ(*,*) ISCR
!
!     Read in the DATAINH input file
!
      READ(LUI) DV
      IF (DV .NE. UV) THEN
        WRITE(6,20) DV,UV
        GO TO 999
      ELSE
        REWIND (UNIT=LUI)
      ENDIF
20    FORMAT(' ERROR:  DATAINH Version No. ',F4.2,' does not match',/,  &
     &       '         UNSAT-H Version No. ',F4.2,'.  Program aborted.')
      READ(LUI) DV,IPLANT,LOWER,NDAYS,NDAY,NPRINT,ITOPBC,ICONVH,MAT,    &
     &          KOPT,KEST,IVAPOR,IDEND,IFDEND,NYEARS,IYS,ISWDIF,ICLOUD, &
     &          NPT,G,MAXPOL,MAXCOE,IEVOPT,NFPET,NSOW,NHRVST,IDTEND,    &
     &          INC,MXROOT,IETOPT,ISHOPT,INMAX,ISTEAD,ILEAP,            &
     &          IRAIN,IHEAT,UPPERH,LOWERH,IFILE,HYFILE,TITLE
      READ(LUI) DMAXBA,DELMAX,DELMIN,RAINIF,RFACT,HIRRI,HDRY,SATK,      &
     &          DRYK,SATC,DRYC,SATTH,DRYTH,AA,B1,B2,TMOIST,LOG10E,      &
     &          DHMAX,DHFACT,QHCTOP,OUTTIM,DMAXHE,HTOP,TGRAD,UP,DOWN,   &
     &          STOPHR,GRAV,TSMEAN,TSAMP,DLAI,BARE,HPR
      READ(LUI) (Z(I),H(I),NTROOT(I),THETA(I),KL(I),C(I),T(I),          &
     &          CHSOIL(I),I=1,NPT),THETAW,THETAD,THETAN,RDF,FPET,       &
     &          MGR,VAPDIF,VC,TSOIL,PETPC,QWLEAK,QHLEAK
      READ(LUI) SWPA,SHPA,SINLAT,COSLAT,TANLAT,TPIY,DAYSEXT,SB,CHW,     &
     &          ALBEDO,ALT,PMB,ZU,ZH,ZM,ZT,D,VK,SEXT,CHA,               &
     &          TCON,CHS,EF,WATDEN,LHV0,CHV,                            &
     &          T0,HSTORE,ISMETH,DHTOL,RHA,IHYS,HYSHPH,SARWA,IPATHA,    &
     &          AIRTOL,HYSTOL,HYSMXH,INDEXB,BBHEAD,HM,THTA,ALFACT
      READ(LUI) (PMFN(I),I=1,NYEARS)
      READ(LUI) (PRFN(I),I=1,NYEARS)
      READ(LUI) PTRANS,PEVAPO
      READ(LUI) NWATER
      IDBEG  = NDAY+1
      DELT   = DELMIN
      NPPT   = NPT-1
      NTOTAL = NINT(24.D0/OUTTIM)
      NSTOP  = NINT(NTOTAL*STOPHR/24.D0)
      DELSUB = 24.D0/NTOTAL
      NGDAYS = NHRVST-NSOW+1
      IF (NSOW .GT. NHRVST) NGDAYS = 365-NSOW+1+NHRVST
      IF (ITOPBC .EQ. 1) H(1) = HTOP
      DO I=1,NPT
        HH(I) = H(I)
      ENDDO
!
!     Calculate internodal distances (ZZ) and node thicknesses (THICK)
!
      ZZ(1) = Z(2)-Z(1)
      THICK(1) = 0.5*ZZ(1)
      THICK(NPT) = 0.5*(Z(NPT)-Z(NPT-1))
      DO I=2,NPT-1
        ZZ(I) = Z(I+1)-Z(I)
        THICK(I) = 0.5*(ZZ(I-1)+ZZ(I))
      ENDDO
      IF (IVAPOR .EQ. 1) THEN
        CALL VOLVAP(1,NPT)
        CALL KVCALC(1,NPT,TSOIL,VC)
      ENDIF
      HSOL(1) = H(1)
      IF (IHYS .GT. 0) THEN
        CALL HYSINI
        CALL HYSSHP(1,NPT,0,KL,C,THETA)
!
!     For restarts, and for new starts on a path other than the primary
!     desorption path, estimate total water storage in the profile
!
        TMOIST = 0.0
        DO I=1,NPT
          TMOIST = TMOIST+(THETAV(I)+THETA(I))*THICK(I)
        ENDDO
        STORE = TMOIST
        TPREV = TMOIST
      ENDIF
!
!     Initialize the precipitation data
!     If IDBEG is greater than IRDAY, go through the data until 
!     the appropriate starting day is reached
!
      IF (IETOPT .EQ. 0) THEN
        IF (NWATER .GT. 0) THEN
          READ(LUI) IRDAY,IRTYPE,EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
          WINDEX = 1
          IF (IRTYPE .EQ. 4) THEN
            NEWPON = RTIME(2)
            IF (NEWPON .EQ. 0) NEWPON = IDEND
            ITPOND = 2
          ENDIF
188       IF (IRDAY .LT. IDBEG .AND. WINDEX .LT. NWATER) THEN
            READ(LUI) IRDAY,IRTYPE,EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
            WINDEX = WINDEX+1
            GO TO 188
          ENDIF
        ENDIF
      ELSE
        METLES = 0
187     READ(LUI) METDAY,TMAX,TMIN,TDEW,SR_MEAS,WIND,CLOUD,IRTYPE,      &
     &      EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
        IF (METDAY .LT. IDBEG) THEN
          TMAXLAST = TMAX
          METLES = 1
          GOTO 187
        ENDIF
        IF (METLES .EQ. 0) TMAXLAST = TMAX
        IF (METDAY .LT. IDTEND) THEN
          READ(LUI) II,RR,TMINNEXT
          BACKSPACE 1
        ELSE
          TMINNEXT = TMIN
        ENDIF
        TAMEAN = (TMAX+TMIN)*.5
        VD_A   = FUNVD(TDEW)
        SVD_A  = FUNVD(TAMEAN)
        RHMEAN = VD_A/SVD_A
        IF (ISHOPT .GT. 0 .AND. IHEAT .EQ. 0) THEN
          IF (ISHOPT .EQ. 1) THEN
            RHA = RHMEAN
          ELSE IF (ISHOPT .EQ. 2) THEN
            TATIME = HOUR+0.5*DELSUB
            CALL AIRTMP(TATIME,TAMEAN,TMAXLAST,TMINNEXT,TA)
            RHA = VD_A/FUNVD(TA)
          ENDIF
          IF (RHA .GT. 1.0) RHA = 1.0
          IF (IEVOPT .EQ. 1) CALL HDRYCALC(RHA)
        ENDIF
        IF (NP .GT. 0) IRDAY = METDAY
      ENDIF
      IF (IHEAT .EQ. 0) THEN
        IF (HDRY .EQ. 0) THEN
          IF (IEVOPT .EQ. 1 .AND. ISHOPT .EQ. 0) CALL HDRYCALC(RHA)
        ENDIF
        RRHA = RHA
        IF (IEVOPT .GT. 1) THEN
          CALL RELHUM(H(1),TS,RHS)
          RRHS = RHS
          PESUB = PEVAPO(NDAY)*FPET(1)
          IF (RRHA .LT. 1) THEN
            EEFLUX = -PESUB*(RRHS-RRHA)/(1-RRHA)
            IF (EEFLUX .GT. 0.0) EEFLUX = 0.0
          ENDIF
        ENDIF
      ENDIF
      HOUR = NTOTAL*OUTTIM
      CALL ZEROA(NP+1,25,RTIME)
      CALL ZEROA(NP+1,25,AMOUNT)
      IF (LOWER .EQ. 2) HSOL(NPT) = H(NPT)
      IF (LOWER .EQ. 4) QLAST = 0.0
      IF (IHEAT .EQ. 1) THEN
        CALL THERMK(1,NPT)
        IF (UPPERH .EQ. 0) THEN
          RLNZT = LOG((ZT-D+ZH)/ZH)
          RLNZU = LOG((ZU-D+ZM)/ZM)
          TATIME = 0.0
          CALL AIRTMP(TATIME,TAMEAN,TMAXLAST,TMINNEXT,TA)
          TTA = TA
          TSMEAN = TS
          CALL BLR(RH)
        ENDIF
      ENDIF
      CALL INTERK(1,NPPT,KEST,UP,DOWN)
!
!     Storing initial values (e.g., H, THETA) as previous values, where
!     a double first letter implies "previous value" (e.g., HH, TTHETA)
!
      CALL RESET(1,NPT,CC,HH,TTHETA,KKLMID,KKVMID,KKVTMID,TT,CCH,       &
     &     KKHMID,TTHETAV,KKLNPT,KKHNPT,KKVTNPT,C,H,THETA,KLMID,KVMID,  &
     &     KVTMID,T,CH,KHMID,THETAV,KL(NPT),KH(NPT),KVT(NPT))
!-----------------------------------------------------------------------
!     Multiple Year Loop  
!        
!     Head values are retained from end of previous year
!     Binary input data file is rewound and repeated each year for
!     steady case, or new PET/MET and precip files are read for unsteady
!     case.
!-----------------------------------------------------------------------
      DO IYEAR=1,NYEARS
      HOUR   = 0.0
      STORE  = TMOIST
      TPREV  = TMOIST
      THPREV = HSTORE
      DHPREV = HSTORE
      IF (IYEAR .GT. 1) NDAY = 0
      IDAY  = NDAY
      IDBEG = NDAY+1
      IF (IYEAR .GT. 1) THEN
        JYEAR = IYS+IYEAR-1
        IF (ISTEAD .EQ. 0) THEN
          ILEAP = 0
          IF (MOD(JYEAR,4) .EQ. 0) THEN
            ILEAP = 1
            IF (JYEAR .GE. 100 .AND. MOD(JYEAR,100) .EQ. 0) ILEAP = 0
            IF (JYEAR .GE. 400 .AND. MOD(JYEAR,400) .EQ. 0) ILEAP = 1 
          ENDIF
        ENDIF
        IF (IYEAR .LT. NYEARS) THEN
          IDEND = 365+ILEAP
          NDAYS = IDEND
        ELSE
          IDEND = IFDEND
          NDAYS = IDTEND
        ENDIF
        IF (ISTEAD .EQ. 0) THEN
          CLOSE (UNIT=LUI)
          OPEN(UNIT=LUX,FILE=PMFN(IYEAR),STATUS='OLD')
          CALL ETBC(0,.FALSE.,LUX)
          CLOSE (UNIT=LUX)
          IF (IPLANT .EQ. 1) CALL PETPART(0,.FALSE.)
          IF (IRAIN .EQ. 0) THEN
            OPEN(UNIT=LUX,FILE=PRFN(IYEAR),STATUS='OLD')
            OPEN(UNIT=LUI,FILE='prec.dat',STATUS='UNKNOWN')
            CLOSE(UNIT=LUI,STATUS='DELETE')
            OPEN(UNIT=LUI,FILE='prec.dat',STATUS='NEW',                 & 
     &        FORM='UNFORMATTED')
          ENDIF
          CALL PRECIPBC(0,.FALSE.,LUX,LUI)
          CLOSE (UNIT=LUX)
          REWIND LUI
          READ(LUI)
        ELSE
          REWIND LUI
          DO I=1,6
            READ(LUI)
          ENDDO
          READ(LUI) PTRANS,PEVAPO
          READ(LUI) NWATER
        ENDIF
        CALL ZEROA(1,25,RTIME)
        CALL ZEROA(1,25,AMOUNT)
        IF (IETOPT .EQ. 0) THEN
          IF (NWATER .GT. 0) THEN
            READ(LUI) IRDAY,IRTYPE,EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
            WINDEX = 1
          ENDIF
        ELSE
          TMAXLAST = TMAX
          READ(LUI) METDAY,TMAX,TMIN,TDEW,SR_MEAS,WIND,CLOUD,IRTYPE,    &
     &        EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
          IF (METDAY .LT. IDTEND) THEN
            READ(LUI) II,RR,TMINNEXT
            BACKSPACE 1
          ELSE
            TMINNEXT = TMIN
          ENDIF
          TAMEAN  = (TMAX+TMIN)*.5
          VD_A   = FUNVD(TDEW)
          SVD_A = FUNVD(TAMEAN)
          RHMEAN = VD_A/SVD_A
          IF (ISHOPT .GT. 0 .AND. IHEAT .EQ. 0) THEN
            IF (ISHOPT .EQ. 1) THEN
              RHA = RHMEAN
            ELSE IF (ISHOPT .EQ. 2) THEN
              TATIME = HOUR+0.5*DELSUB
              CALL AIRTMP(TATIME,TAMEAN,TMAXLAST,TMINNEXT,TA)
              RHA = VD_A/FUNVD(TA)
            ENDIF
            IF (RHA .GT. 1.0) RHA = 1.0
            RRHA = RHA
            IF (IEVOPT .EQ. 1) CALL HDRYCALC(RHA)
            IF (IEVOPT .GT. 1) THEN
              CALL RELHUM(H(1),TS,RHS)
              RRHS = RHS
            ENDIF
          ENDIF
          IF (NP .GT. 0) IRDAY = METDAY
        ENDIF
      ENDIF
!
!     Write initial conditions to binary output file
!
      IF (NYEARS .EQ. 1) THEN
        RFILE(NCHR:NCHR+3) = '.res'
      ELSE
        RFILE(NCHR:NCHR+3) = '0000'
        NYEAR = IYS+IYEAR-1
        WRITE(UNIT=AYEAR,FMT='(I)') NYEAR
        CALL CLIP(AYEAR)
        NCHRY = INDEX(AYEAR," ")-1
        RFILE(NCHR+4-NCHRY:NCHR+3) = AYEAR(1:NCHRY)
        RFILE(NCHR+4:NCHR+7) = '.res'
      ENDIF
      IRL = IRLMOD*MAX(110,(9*IHEAT+IPLANT+5)*NPT+45)
      OPEN (UNIT=LUB,FILE=RFILE,STATUS='UNKNOWN')
      CLOSE(UNIT=LUB,STATUS='DELETE')
      OPEN (UNIT=LUB,FILE=RFILE,STATUS='NEW',FORM='UNFORMATTED',        &
     &  ACCESS='DIRECT',RECL=IRL)
      WRITE(LUB,REC=1) NPT,IPLANT,IHEAT,(Z(I),I=1,NPT),UV,STOPHR,TSOIL, &
     &  IFILE,SDATE,STIME,TITLE,IDEND,NPRINT,IEVOPT,NTOTAL,IDBEG,IVAPOR,&
     &  ICLOUD,IETOPT,ISHOPT,UPPER,UPPERH,LOWER,LOWERH,IYEAR
!
!     Zero daily cumulative totals and write initial conditions
!     to binary output file
!
      HOUR = 24.0
      CALL ZEROI(DAYAST,DAYSTP,DAYSDH,DAYUBC)
      CALL ZEROR(DAYE,DAYINF,DAYIRR,DAYPE,DAYPT,DAYRAN,DAYRUN,DAYTIM,   &
     &  DAYTRA,DAYHBE,DAYLE,DAYRN,DAYSHF,DAYSEN)
      DAYQHW0 = 0.0
      CALL ZEROA(1,NPT,DAYSNK)
      CALL ZEROA(1,NPT,DAYQL)
      CALL ZEROA(1,NPT,DAYQVH)
      CALL ZEROA(1,NPT,DAYQVT)
      CALL ZEROA(1,NPT,DAYQHC)
      CALL ZEROA(1,NPT,DAYQHV)
      CALL ZEROA(1,NPT,DAYQHW)
      IF (IPLANT .EQ. 0) THEN
        IF (IHEAT .EQ. 0) THEN
          WRITE(LUB,REC=2) IDAY,HOUR,(H(J),THETA(J),DAYQL(J),DAYQVH(J), &
     &    J=1,NPT),DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,                    &
     &    TPREV,TMOIST,DAYSTP,DAYPE,DAYPT,DAYTIM,DAYAST,DAYUBC,RHMEAN,  &
     &    TAMEAN,HDRY,DAYRN,DAYSHF,DAYSEN,DHPREV,HSTORE,DAYQHW0,        &
     &    TA,TMAX,TMIN,VD_A,WIND,CLOUD,SR_MEAS,DAYLE,DAYSDH,DAYHBE
        ELSE
          WRITE(LUB,REC=2) IDAY,HOUR,(H(J),THETA(J),DAYQL(J),DAYQVH(J), &
     &    DAYQVT(J),T(J),J=1,NPT),DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,     &
     &    TPREV,TMOIST,DAYSTP,DAYPE,DAYPT,DAYTIM,DAYAST,DAYUBC,RHMEAN,  &
     &    TAMEAN,HDRY,DAYRN,DAYSHF,DAYSEN,DHPREV,HSTORE,DAYQHW0,        &
     &    TA,TMAX,TMIN,VD_A,WIND,CLOUD,SR_MEAS,DAYLE,DAYSDH,DAYHBE,     &
     &    (DAYQHC(J),DAYQHW(J),DAYQHV(J),J=1,NPT)
        ENDIF
      ELSE
        WRITE(LUB,REC=2) IDAY,HOUR,(H(J),THETA(J),DAYQL(J),DAYQVH(J),   &
     &  DAYSNK(J),J=1,NPT),DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,            &
     &  TPREV,TMOIST,DAYSTP,DAYPE,DAYPT,DAYTIM,DAYAST,DAYUBC,RHMEAN,    &
     &  TAMEAN,HDRY,DAYSDH
      ENDIF
      IREC = 2
!
!     Initialize yearly cumulatives
!
      APLIED = 0.0
      OLDRAT = 0.0
      TETRAN = 0.0
      TPET   = 0.0
      CALL ZEROI(TAST,TSTP,TSDH,TUBC)
      CALL ZEROR(TE,TINF,TIRR,TPE,TPT,TRAN,TRUN,TTIM,                   &
     &  TTRA,xxx,TLE,TRN,TSHF,TSEN)
      TQHW0 = 0.0
      CALL ZEROA(1,NPT,TSNK)
      CALL ZEROA(1,NPT,TQL)
      CALL ZEROA(1,NPT,TQVH)
      CALL ZEROA(1,NPT,TQVT)
      CALL ZEROA(1,NPT,TQHC)
      CALL ZEROA(1,NPT,TQHV)
      CALL ZEROA(1,NPT,TQHW)
!-----------------------------------------------------------------------
!     Day Loop
!
!     IDAY   = Day number
!     IDBEG = Beginning day
!     IDEND = Ending day
!-----------------------------------------------------------------------
      DO IDAY=IDBEG,IDEND
      TPREV  = TMOIST
      DHPREV = HSTORE
      IF (LOWER .EQ. 3) QLAST = QWLEAK(IDAY)/24.
      TOTALR = 0.
      IF (IRTYPE .NE. 4) THEN
        IRRI  = 0
        IPREC = 0
      ENDIF
      CALL ZEROA(1,24,RATE)
!
!     Check for rain on current day and set times
!
      IF (IETOPT .EQ. 0) THEN
        IF (IRDAY .NE. IDAY) GO TO 330
      ELSE
        IF (NP .EQ. 0) GO TO 330
      ENDIF
      IF (IRTYPE .EQ. 1) THEN
        IPREC = 1
      ELSE
        IRRI = 1
      ENDIF
!
!     Determines the water input rate (cm/h)
!     If IRTYPE = 1, input is PPT, if IRTYPE = 2, input is IRRI
!
      IF (IRTYPE .EQ. 1 .OR. IRTYPE .EQ. 2) THEN
        IF (NP .EQ. 0) GO TO 330
        DO I=1,NP-1
          TOTALR  = TOTALR+AMOUNT(I)
        ENDDO
        IF (IRTYPE .EQ. 2) TOTALR = TOTALR*EFICEN
        APLIED = APLIED+TOTALR
        J = 1
        DO I=1,24
          IF ((I-1) .GE. RTIME(1) .AND. (I-1) .LT. RTIME(NP)) THEN
            IF ((I-1) .EQ. RTIME(J)) THEN
              RATE(I) = AMOUNT(J)/(RTIME(J+1)-RTIME(J))
              J = J+1
            ELSE
              RATE(I) = RATE(I-1)
            ENDIF
          ENDIF
        ENDDO
      ELSE IF (IRTYPE .EQ. 3) THEN 
        HPOND = AMOUNT(1)
      ELSE IF (IRTYPE .EQ. 4) THEN
        HPOND = AMOUNT(ITPOND)
        IHPOND = 1
        IF (IDAY .GT. NEWPON) THEN
          ITPOND = ITPOND+1
          NEWPON = RTIME(ITPOND)
          IF (NEWPON .EQ. 0) NEWPON = IDEND
        ENDIF
      ENDIF
330   IF (IHEAT .EQ. 0) THEN
!
!     If IPLANT = 0, PTRANS = 0 and PET becomes PEVAPO
!     If IPLANT = 1, PET becomes PTRANS + PEVAPO
!
        PE = PEVAPO(IDAY)
        PT = PTRANS(IDAY)
      ELSE
        IF (UPPERH .EQ. 0) CALL NETRAD(1,1)
      ENDIF
      IF (IPLANT .GT. 0) THEN
        IF (IDAY .GE. NSOW) THEN
          NGROW = IDAY-NSOW+1
        ELSE
          NGROW = 365+IDAY-NSOW+1
        ENDIF
!
!     Determine rooting depth and density
!
        IF (NGROW .LE. NGDAYS) CALL ROOTS(NGROW,MDEPTH)
      ENDIF
!
!     Initialize daily cumulative totals
!
      HOUR   = 0.0
      PTSUB  = 0.0
      PESUB  = 0.0
      TSTEPB = 0.0
      CALL ZEROI(DAYAST,DAYSTP,DAYSDH,DAYUBC)
      CALL ZEROR(DAYE,DAYINF,DAYIRR,DAYPE,DAYPT,DAYRAN,DAYRUN,DAYTIM,   &
     &  DAYTRA,DAYHBE,DAYLE,DAYRN,DAYSHF,DAYSEN)
      DAYQHW0 = 0.0
      CALL ZEROA(1,NPT,DAYSNK)
      CALL ZEROA(1,NPT,DAYQL)
      CALL ZEROA(1,NPT,DAYQVH)
      CALL ZEROA(1,NPT,DAYQVT)
      CALL ZEROA(1,NPT,DAYQHC)
      CALL ZEROA(1,NPT,DAYQHV)
      CALL ZEROA(1,NPT,DAYQHW)
!-----------------------------------------------------------------------
!     NTOTAL Loop
!
!     Program will remain in this loop until 1 day has been simulated.
!-----------------------------------------------------------------------
      DO N=1,NTOTAL
        SPREV  = TMOIST
        SHPREV = HSTORE
        NHR    = INT(1.0+(N-1)*24.0/NTOTAL)
        IF (IPLANT .EQ. 1) PTSUB = PT*FPET(NHR)
        IF (IEVOPT .GT. 0 .AND. IHEAT .EQ. 0) PESUB = PE*FPET(NHR)
        IF (IPREC .EQ. 1) PRAIN = RATE(NHR)
        IF (IRRI .EQ. 1) PIRRI = RATE(NHR)
        IF (RATE(NHR) .EQ. 0.) THEN
          ISATUK = 0
        ELSE
          ISATUK = 1
          IF (OLDRAT .GT. 0.) ISATUK = 2
        ENDIF
        IF (ISHOPT .EQ. 2 .AND. IHEAT .EQ. 0) THEN
          TATIME = HOUR+0.5*DELSUB
          CALL AIRTMP(TATIME,TAMEAN,TMAXLAST,TMINNEXT,TA)
          RHA = VD_A/FUNVD(TA)
          IF (RHA .GT. 1.0) RHA = 1.0
          IF (IEVOPT .EQ. 1) CALL HDRYCALC(RHA)
        ENDIF 
!
!     Initialize total values for each OUTTIM (output time)
!
      CALL ZEROI(SUBAST,SUBSTP,SUBSDH,SUBUBC)
      CALL ZEROR(SUBE,SUBINF,SUBIRR,SUBPE,SUBPT,SUBRAN,SUBRUN,SUBTIM,   &
     &  SUBTRA,SUBHBE,SUBLE,SUBRN,SUBSHF,SUBSEN)
      SUBQHW0 = 0.0
      CALL ZEROA(1,NPT,SUBSNK)
      CALL ZEROA(1,NPT,SUBQL)
      CALL ZEROA(1,NPT,SUBQVH)
      CALL ZEROA(1,NPT,SUBQVT)
      CALL ZEROA(1,NPT,SUBQHC)
      CALL ZEROA(1,NPT,SUBQHV)
      CALL ZEROA(1,NPT,SUBQHW)
!-----------------------------------------------------------------------
!     Start of the DELSUB loop. The program remains here until DELSUB
!     time has elapsed. (i.e. TSUB = DELSUB)
!     TSTEPB = time of day (h) at beginning of time step
!     TSTEPM = time of day (h) at midpoint of time step
!     TSTEPE = time of day (h) at end of time step
!-----------------------------------------------------------------------
      TSUB = 0.0
32649 IF (DELSUB-TSUB .LE. TIMRES) GO TO 32648
      IF (DELSUB-TSUB .LT. DELMIN) DELT = DELSUB-TSUB
      DELSAV = DELT
      TSTEPB = HOUR+TSUB
      TSUB   = TSUB+DELSAV
      CALL INTERK(1,NPPT,KEST,UP,DOWN)
!
!     When UPPER = 0, Water flux surface bc
!     When UPPER = 1, constant head surface bc
!
      UPPER = 0
      NSINK = 0
      SFLUX = 0.0
      IF (IPLANT .EQ. 1) THEN
        IF (NGROW .GT. 0 .AND. NGROW .LE. NGDAYS) NSINK = 1
      ENDIF
      IF (IHEAT .EQ. 1 .AND. UPPERH .EQ. 0) THEN
        SSVD_S = FUNVD(TT(1))
        CALL RELHUM(HH(1),TT(1),RHUMID)
        VVD_S = SSVD_S*RHUMID
      ENDIF
      IF (IRRI .EQ. 1 .AND. IRTYPE .GT. 2) THEN
        IF (HPOND .GT. 0) THEN
          H(1) = -HPOND
          UPPER = 1
        ELSE
          IHPOND = 0
          IRRI   = 0
        ENDIF
      ELSE
!
!     During watering, surface evaporation set to zero
!
        IF (ISATUK .GT. 0) THEN
          IF (ISATUK .EQ. 1) THEN
!
!     For the first time step of a rainfall event,
!     DELT is reduced by RAINIF
!
            TSUB  = TSUB-DELT
            DELT  = DELT*RAINIF
            IF (DELT .LT. DELMIN) DELT = DELMIN
            TSUB = TSUB+DELT
          ENDIF
          NSINK = 0
          SFLUX = RATE(NHR)
        ENDIF
      ENDIF
!
!     When no infiltration, surface flux set to the evaporation rate
!
      IF (IEVOPT .EQ. 1 .AND. IHEAT .EQ. 0) THEN
        IF (ISATUK .EQ. 0 .AND. H(1) .LE. HDRY) SFLUX = -PESUB
      ENDIF
!
!     Calculate plant sink term
!
      IF (IPLANT .GT. 0) CALL TRSINK(MDEPTH,PTSUB,NSINK,SINK,TRANSP)
      ISAVUP = UPPER
      IF (ITOPBC .EQ. 1) THEN
        UPPER  = 1
        ISAVUP = UPPER
      ENDIF
!-----------------------------------------------------------------------
!     "Decrease DELT" loop when solution unacceptable
!-----------------------------------------------------------------------
150   DELSAV = DELT
      IN = 0
      TSTEPM = TSTEPB+0.5*DELSAV
      TSTEPE = TSTEPB+DELSAV
      IF (ISHOPT .EQ. 3 .OR. IHEAT .EQ. 1) THEN
        CALL AIRTMP(TSTEPE,TAMEAN,TMAXLAST,TMINNEXT,TA)
        IF (ISHOPT .EQ. 3) THEN
          RHA = VD_A/FUNVD(TA)
          IF (RHA .GT. 1.0) RHA = 1.0
        ENDIF
      ENDIF
!
!     When LOWER = 5, specify bottom boundary head as function of time
!
      IF (LOWER .EQ. 5 .AND. INDEXB .GT. 0) THEN
        POTTIM = IDAY-1+(TSTEPE)/24.0
        IB = INDEXB
        DO I=2,INDEXB
          IF (POTTIM .GT. BBHEAD(I-1,1) .AND.                           &
     &        POTTIM .LE. BBHEAD(I,1)) IB = I
        ENDDO
        BBTIME = (POTTIM-BBHEAD(IB-1,1))/(BBHEAD(IB,1)-BBHEAD(IB-1,1))
        H(NPT) = BBTIME*BBHEAD(IB,2)+(1-BBTIME)*BBHEAD(IB-1,2)
        HSOL(NPT) = H(NPT)
      ENDIF
      ITSUCC = 0
!-----------------------------------------------------------------------
!     Predictor-Corrector Loop (while IN .LT. INMAX )
!
!     This implicit scheme solves the partial differential equation
!     for the change in head
!-----------------------------------------------------------------------
32626 IF (IN .GE. INMAX .OR. ITSUCC .EQ. 1) GO TO 32625
      IN = IN+1
      DEL = 1.0/DELSAV
      IF (IHEAT .EQ. 1 .AND. UPPERH .EQ. 0) CALL NETRAD(2,IN)
      IF (IN .GT. 1) THEN
        IF (IVAPOR .EQ. 1) THEN
          CALL VOLVAP(1,NPT)
          CALL KVCALC(1,NPT,TSOIL,VC)
        ENDIF
        IF (IHEAT  .EQ. 1) CALL THERMK(1,NPT)
        CALL INTERK(1,NPPT,KEST,UP,DOWN)
      ENDIF
!
!     Calculate coefficients for simulation of change in head
!
556   DO II=2,NPPT
        I      = II-1
        III    = II+1
        Z31    = Z(III)-Z(I)
        FZZ1   = 1./(Z31*ZZ(I))
        FZZ2   = 1./(Z31*ZZ(II))
        ZK1    = (KLMID(I)+KVMID(I))*FZZ1
        ZK2    = (KLMID(II)+KVMID(II))*FZZ2
        A1(II) = ZK1
        A3(II) = ZK2
        PG1    = (HH(II)-HH(I))/ZZ(I)
        PG2    = (HH(III)-HH(II))/ZZ(II)
        IF (ISMETH .EQ. 0) THEN
          AVC = 0.5*(C(II)+CC(II))
          IF (IHYS .GT. 0) THEN
!
!       when the path reverses (JH .NE. JJH), use the forward
!       capacity C for AVC for that time step
!
            IF (JH(II) .NE. JJH(II)) AVC = C(II)
          ENDIF
          A2(II) = AVC*DEL-ZK1-ZK2
          ZY(II) = AVC*HH(II)*DEL-SINK(II)-(KKLMID(II)*(PG2+G)          &
     &    +KKVMID(II)*PG2-KKLMID(I)*(PG1+G)-KKVMID(I)*PG1+G*(KLMID(II)  &
     &    -KLMID(I)))/Z31
        ELSE
          A2(II) = C(II)*DEL-ZK1-ZK2
          ZY(II) = (TTHETA(II)-THETA(II)+C(II)*H(II))*DEL               &
     &    -SINK(II)-(KKLMID(II)*(PG2+G)+KKVMID(II)                      &
     &    *PG2-KKLMID(I)*(PG1+G)-KKVMID(I)*PG1+G*(KLMID(II)             &
     &    -KLMID(I)))/Z31
        ENDIF
      ENDDO
      IF (IHEAT .EQ. 1 .AND. IVAPOR .EQ. 1) THEN
        DO II=2,NPPT
          I      = II-1
          III    = II+1
          Z31    = Z(III)-Z(I)
          ZY(II) = ZY(II)-((-KVTMID(II)*(T(III)-T(II))-KKVTMID(II)*     &
     &             (TT(III)-TT(II)))/ZZ(II)+(KVTMID(I)*(T(II)-T(I))+    &
     &             KKVTMID(I)*(TT(II)-TT(I)))/ZZ(I))/Z31
        ENDDO
      ENDIF
!-----------------------------------------------------------------------
!     Surface node calculations
!-----------------------------------------------------------------------
      IF (IHEAT .EQ. 0) THEN
        IF (IEVOPT .GT. 1 .AND. ISATUK .EQ. 0) THEN
          IF (IEVOPT .EQ. 2) THEN
            IF (IN .EQ. 1) CALL RELHUM(HH(1),TS,RHS)
          ELSE
            CALL RELHUM(H(1),TS,RHS)
          ENDIF
          EFLUX = 0.0
          IF (RHA .LT. 1.0) EFLUX = -PESUB*(RHS-RHA)/(1-RHA)
          IF (EFLUX .GT. 0.) EFLUX = 0.0
          SFLUX = 0.5*(EEFLUX+EFLUX)
          IF (SFLUX .GT. 0.0) SFLUX = 0.0
          IF (IEVOPT .EQ. 3) THEN
            DEDH = 0.0
            CAMLHS = 0.0
            IF (EFLUX .LT. 0.0) THEN
              DEDH   = 0.5*PESUB*MGR*RHS/(TS*(1-RHA))
              SFLUX  = SFLUX-DEDH*H(1)
              CAMLHS = DEDH
            ENDIF
          ENDIF
        ENDIF
      ELSE
        VDSLHS = 0.0
        VDSRHS = 0.0
        IF (ISATUK .EQ. 0 .AND. UPPERH .EQ. 0) THEN
          SVD_S = FUNVD(TS)
          CALL RELHUM(H(1),TS,RHUMID)
          VD_S = SVD_S*RHUMID
!
!     Because water and vapor density are in g/cm3 and RH is in s/m, 
!     the calculated evaporation rate is multiplied by 3.6E5 cm s/m hr 
!     to yield a SFLUX value with units of cm/hr.
!
!     The VDSLHS term is used to modify the A2 term in the coefficient
!     matrix to account for the x-terms in the expansion of 
!     exp(x) to 1+x+x^2/2&... in the calculation of VD_S.
!
          TSMEAN = 0.5*(TT(1)+TS)
          CALL BLR(RH)
          TERMH1 = -(H(1)+OSMPOT)*MGR/TS
          TERM   = TERMH1
          VDSRHS = 1.0-OSMPOT*MGR/TS
          NX     = 4+3E-6*H(1)
          DO I = 2,NX
            TERM   = TERM*TERMH1/I
            VDSRHS = VDSRHS+TERM
          END DO
          VDSRHS = SVD_S*VDSRHS
          SFLUX  = 3.6D5*(VD_A-((VVD_S+VDSRHS)/2.0))/(WATDEN*RH)
          VDSLHS = (MGR/TS)*3.6D5*SVD_S/(2*WATDEN*RH)
!
!     Positive water fluxes (to the surface) caused by vapor density
!     gradients are not allowed to occur (i.e., no dew additions)
!
          IF (VD_A .GE. 0.5*(VVD_S+VD_S)) THEN
            SFLUX  = 0.0
            VDSLHS = 0.0
          ENDIF
        ENDIF
      ENDIF
      IF (UPPER .EQ. 1) THEN
!
!     Surface node equation for constant head b.c.
!
        ISTART = 2
        ZY(2)  = ZY(2)-A1(2)*H(1)
        A1(2)  = 0.0D0
      ELSE
!
!     Surface node equation for flux b.c.
!
        Z1     = ZZ(1)
        Z11    = (KLMID(1)+KVMID(1))/(Z1*Z1)
        A1(1)  = 0.0
        A3(1)  = Z11
        PG2    = (HH(2)-HH(1))/Z1
        IF (ISMETH .EQ. 0) THEN
          AVC = (C(1)+CC(1))*0.5
          IF (IHYS .GT. 0) THEN
            IF (JH(1) .NE. JJH(1)) AVC = C(1)
          ENDIF
          A2(1) = AVC*DEL-Z11-2*VDSLHS/Z1-2*CAMLHS/Z1
          ZY(1) = AVC*DEL*HH(1)-(KKLMID(1)*(PG2+G)+KKVMID(1)            &
     &            *PG2+G*KLMID(1)-2.0*SFLUX)/Z1
        ELSE
          A2(1) = C(1)*DEL-Z11-2*VDSLHS/Z1-2*CAMLHS/Z1
          ZY(1) = (TTHETA(1)-THETA(1)+C(1)*H(1))*DEL                    &
     &            -(KKLMID(1)*(PG2+G)+KKVMID(1)                         &
     &            *PG2+G*KLMID(1)-2.0*SFLUX)/Z1
        ENDIF
        IF (IHEAT .EQ. 1 .AND. IVAPOR .EQ. 1) ZY(1) = ZY(1)-            &
     &    ((-KVTMID(1)*(T(2)-T(1))-KKVTMID(1)*(TT(2)-TT(1)))/Z1)/Z1
        ISTART = 1
      ENDIF
!
!     LOWER: 1=unit grad.; 2=cons. head; 3=specified flux; 4=imperm.;
!            5=variable head specified by user
!
      IF (LOWER .NE. 2 .AND. LOWER .NE. 5) THEN
        ZN      = ZZ(NPPT)
        ZK1     = (KLMID(NPPT)+KVMID(NPPT))/(ZN*ZN)
        A3(NPT) = 0.0
        IF (LOWER .EQ. 1) QLAST = 0.5*(KL(NPT)+KKLNPT)
        A1(NPT) = ZK1
        PG1     = (HH(NPT)-HH(NPPT))/ZN
        IF (ISMETH .EQ. 0) THEN
          AVC = 0.5*(C(NPT)+CC(NPT))
          IF (IHYS .GT. 0) THEN
            IF (JH(NPT) .NE. JJH(NPT)) AVC = C(NPT)
          ENDIF
          A2(NPT) = AVC*DEL-ZK1
          ZY(NPT) = AVC*HH(NPT)*DEL-SINK(NPT)                           &
     &             -(2.0*QLAST-KKLMID(NPPT)*(PG1+G)                     &
     &             -KKVMID(NPPT)*PG1-G*KLMID(NPPT))/ZN
        ELSE
          A2(NPT) = C(NPT)*DEL-ZK1
          ZY(NPT) = (TTHETA(NPT)-THETA(NPT)+C(NPT)*H(NPT))*DEL-SINK(NPT)&
     &              -(2.0*QLAST-KKLMID(NPPT)*(PG1+G)                    &
     &              -KKVMID(NPPT)*PG1-G*KLMID(NPPT))/ZN
        ENDIF
        IF (IHEAT .EQ. 1 .AND. IVAPOR .EQ. 1) THEN
          ZY(NPT) = ZY(NPT)-((-KVT(NPT)-KKVTNPT)*TGRAD                  &
     &              +(KVTMID(NPPT)*(T(NPT)-T(NPPT))                     &
     &              +KKVTMID(NPPT)*(TT(NPT)-TT(NPPT)))/ZN)/ZN
        ENDIF
        NPPT = NPT
      ELSE
        ZY(NPPT) = ZY(NPPT)-A3(NPPT)*H(NPT)
        A3(NPPT) = 0.0D0
      ENDIF
!
!     Solution of tridiagonal matrix
!
      CALL TRIDAG(ISTART,NPPT,A1,A2,A3,ZY,HSOL)
!-----------------------------------------------------------------------
!     Solution estimate complete for this iteration
!-----------------------------------------------------------------------
      NPPT = NPT-1
!
!     Reduce time step if predicted head exceeds head lower bound (HLB)
!
      IGO150 = 0
      DO I=1,NPT
        IF (HSOL(I) .LT. HLB) THEN
          IF (DELT .LE. DELMIN) THEN
            WRITE(*,*) ' Head exceeded lower bound using DELMIN'
            WRITE(*,*) ' Simulation stopped in UNSATH'
            GOTO 999
          ENDIF
          TSUB = TSUB-DELT
          DELT = 0.5*DELT
          IF (DELT .LT. DELMIN) DELT = DELMIN
          TSUB = TSUB+DELT
          IGO150 = 1
          GO TO 795
        ENDIF
      ENDDO
!
!     Reduce time step if change in head of any node exceeds DHMAX
!
      IF (DHMAX .GT. 0) THEN
        IF (DELT .GT. DELMIN .AND. ISATUK .NE. 1) THEN
          DO I=1,NPT
            IF (ABS(HH(I)-HSOL(I)) .GE. DHMAX) THEN
              TSUB = TSUB-DELT
              DELT = DELT*DHFACT
              IF (DELT .LT. DELMIN) DELT = DELMIN
              TSUB = TSUB+DELT
              SUBSDH = SUBSDH+1
              IGO150 = 1
              GO TO 795
            ENDIF
          ENDDO
        ENDIF
      ENDIF
!
!     Check whether surface b.c. changed due to solution with given
!     flux.  If so, repeat as constant surface head b.c. 
!
      IF (UPPER .EQ. 0) THEN
        IF (ISATUK .GT. 0) THEN
          IF (HSOL(1) .LT. HIRRI) THEN
            H(1)     = HIRRI
            HSOL(1)  = HIRRI
            UPPER    = 1
            IF (IHYS .GT. 0) CALL HYSSHP(1,1,0,KL,C,THETA)
          ENDIF
        ELSE IF (IHEAT .EQ. 0 .AND. IEVOPT .EQ. 1) THEN
          IF (SFLUX .LT. 0.0 .AND. HSOL(1) .GT. HDRY) THEN
            H(1)     = HDRY
            HSOL(1)  = HDRY
            KL(1)    = DRYK
            C(1)     = DRYC
            THETA(1) = DRYTH
            UPPER    = 1
          ENDIF
        ENDIF 
        IF (UPPER .EQ. 1) THEN
          IF (IN .GT. 1) THEN
            CALL RESET(2,NPT,C,H,THETA,KLMID,KVMID,KVTMID,T,CH,KHMID,   &
     &      THETAV,KL(NPT),KH(NPT),KVT(NPT),CC,HH,TTHETA,KKLMID,KKVMID, &
     &      KKVTMID,TT,CCH,KKHMID,TTHETAV,KKLNPT,KKHNPT,KKVTNPT)
            IN = 1
          ENDIF
          IF (IVAPOR .EQ. 1) THEN
            CALL KVCALC(1,1,TSOIL,VC)
            CALL VOLVAP(1,1)
          ENDIF
          IF (IHEAT .EQ. 1) CALL THERMK(1,1)
          CALL INTERK(1,1,KEST,UP,DOWN)
          GO TO 556
        ENDIF
      ENDIF
!
!     Passed check on upper boundary condition
!
!     Test whether head changes during this iteration are less than
!     DHTOL.  If successful, ITSUCC is set to 1 and the latest head
!     solution is acceptable
!
      IF (DHTOL .GT. 0) THEN
        ITSUCC = 0
        IF (IN .GT. 1 .AND. IN .LT. INMAX) THEN
          HFDH = 0.0
          DO I=1,NPT
            IF (H(I) .NE. 0) THEN
              FDH = ABS( (HSOL(I)-H(I)) / H(I) )
              IF (HFDH .LT. FDH) HFDH = FDH
            ENDIF
          END DO
          IF (HFDH .LT. DHTOL) ITSUCC = 1
        ENDIF
      ENDIF
      DO I=1,NPT
        H(I) = HSOL(I)
      ENDDO
      IF (IHYS .EQ. 0) THEN
        IF (IN .LT. INMAX) THEN
         CALL POLYKH(1,NPT,H,KL)
         CALL RETENT(1,NPT,H,1,THETA,1,C)
        ELSE
          CALL RETENT(1,NPT,H,1,THETA,0,C)
        ENDIF
      ELSE
        CALL HYSSHP(1,NPT,0,KL,C,THETA)
      ENDIF
      IF (IVAPOR .EQ. 1) CALL VOLVAP(1,NPT)
      IF (IHEAT .EQ. 1) THEN
        IGO150 = 0
        CALL FLUX(RH)
        IF (IN .LT. INMAX) CALL KVCALC(1,NPT,TSOIL,VC)
        CALL THERMK(1,NPT)
        CALL HEATFLOW(RH,1)
        IF (IGO150 .EQ. 1) GO TO 795
        IF (IVAPOR .EQ. 1) CALL VOLVAP(1,NPT)
        CALL THERMK(1,NPT)
      ENDIF
      GO TO 32626
!-----------------------------------------------------------------------
!     End of the Predictor-Corrector Loop 
!-----------------------------------------------------------------------
32625 SUBAST = SUBAST+1
      IGO150 = 0
!
!     Fluxes are calculated at all column segments
!
      CALL FLUX(RH)
      IF (IHEAT .EQ. 1) CALL HEATFLOW(RH,2)
!
!     Check for convergence on THETA
!
      IF (ISWDIF .EQ. 0) THEN 
        DIFMAX = 0.0
        DO I=1,NPT
          THDIF = ABS((THETA(I)-TTHETA(I))/TTHETA(I))
          DIFMAX = MAX(DIFMAX,THDIF)
        ENDDO
        CALL DELCHK
        IF (IGO150 .EQ. 1) GO TO 795
      ENDIF
!
!     Estimate total water storage in the profile (SMOIST) and the
!     mass balance error (DIFMAX)
!
      SMOIST = 0.0
      DO I=1,NPT
        SMOIST = SMOIST+(THETA(I)+THETAV(I))*THICK(I)
      ENDDO
      DIFMAX = TMOIST+SSFLUX-TRANSP*DELSAV-QT(NPT)-SMOIST
!
!     Reduce the time step if the mass balance error is not within 
!     limits (when ISWDIF .EQ. 1)
!
      IF (ISWDIF .EQ. 1) CALL DELCHK
795   IF (IGO150 .EQ. 1) THEN
        CALL RESET(1,NPT,C,H,THETA,KLMID,KVMID,KVTMID,T,CH,KHMID,       &
     &    THETAV,KL(NPT),KH(NPT),KVT(NPT),CC,HH,TTHETA,KKLMID,KKVMID,   &
     &    KKVTMID,TT,CCH,KKHMID,TTHETAV,KKLNPT,KKHNPT,KKVTNPT)
        UPPER = ISAVUP
        GO TO 150
      ENDIF
!-----------------------------------------------------------------------
!     End of convergence loop for water flow
!     Time step successful, all conditions satisfied.
!     Old water flow parameters are reset to new values.
!     DELSUB water and heat totals are updated.
!-----------------------------------------------------------------------
      IF (ISATUK .EQ. 1) ISATUK = 2
      SUBSTP = SUBSTP+1
      IF (ISAVUP .EQ. 0 .AND. UPPER .EQ. 1) SUBUBC = SUBUBC+1
      SUBTIM = SUBTIM+DELSAV
      DO I=1,NPT
        SUBQL(I)  = SUBQL(I)+QL(I)
      ENDDO
      IF (IVAPOR .EQ. 1) THEN
        DO I=1,NPT
          SUBQVH(I) = SUBQVH(I)+QVH(I)
          SUBQVT(I) = SUBQVT(I)+QVT(I)
        ENDDO
      ENDIF
      SUBTRA = SUBTRA+TRANSP*DELSAV
      DO I=1,MDEPTH
        SUBSNK(I) = SUBSNK(I)+SINK(I)*DELSAV*THICK(I)
      ENDDO
      IF (SSFLUX .GT. 0.0) SUBINF = SUBINF+SSFLUX
      IF (SSFLUX .LT. 0.0) SUBE = SUBE-SSFLUX
      SUBPE = SUBPE+PESUB*DELSAV
      SUBPT = SUBPT+PTSUB*DELSAV
      IF (RATE(NHR) .GT. 0.0 .OR. HPOND .GT. 0.0) THEN
        SUBRAN = SUBRAN+DELSAV*PRAIN
        SUBIRR = SUBIRR+DELSAV*PIRRI
        RUNOFF = 0.0
        IF (PIRRI*DELSAV .GT. SSFLUX) RUNOFF = PIRRI*DELSAV-SSFLUX
        IF (PRAIN*DELSAV .GT. SSFLUX) RUNOFF = PRAIN*DELSAV-SSFLUX
        IF (RUNOFF .GT. 0.0) SUBRUN = SUBRUN+RUNOFF
      ENDIF
      KKLNPT  = KL(NPT)
      KKVTNPT = KVT(NPT)
      IF (IHYS .EQ. 0) THEN
        CALL RETENT(1,NPT,H,0,THETA,1,C)
        CALL POLYKH(1,NPT,H,KL)
      ELSE
        CALL HYSSHP(1,NPT,1,KL,C,THETA)
      ENDIF
!-----------------------------------------------------------------------
!     If IHEAT = 1, the HEATFLOW subroutine is called to update 
!     the DELSUB total values for heat-related terms
!-----------------------------------------------------------------------
      IF (IHEAT .EQ. 1) THEN
        CALL HEATFLOW(RH,3)
        HSTORE = SHSTORE
        KKHNPT = KH(NPT)
        CALL THERMK(1,NPT)
      ENDIF
      IF (IVAPOR .EQ. 1) THEN
        CALL VOLVAP(1,NPT)
        CALL KVCALC(1,NPT,TSOIL,VC)
      ENDIF
!
!     Store final values (e.g., H) as previous values (e.g., HH)
!
      CALL RESET(1,NPT,CC,HH,TTHETA,KKLMID,KKVMID,KKVTMID,TT,CCH,       &
     &     KKHMID,TTHETAV,KKLNPT,KKHNPT,KKVTNPT,C,H,THETA,KLMID,KVMID,  &
     &     KVTMID,T,CH,KHMID,THETAV,KL(NPT),KH(NPT),KVT(NPT))
      EEFLUX = EFLUX
      RRHS   = RHS
      SSFLUX = 0.0
      TMOIST = SMOIST
      TTA    = TA
      GO TO 32649
!----------------------------------------------------------------------
!     End of the DELSUB Loop.  (TSUB = DELSUB)
!     Daily totals are updated
!----------------------------------------------------------------------
32648 DAYAST = DAYAST+SUBAST
      DAYE   = DAYE+SUBE
      DAYINF = DAYINF+SUBINF
      DAYIRR = DAYIRR+SUBIRR
      DAYPE  = DAYPE+SUBPE
      DAYPT  = DAYPT+SUBPT
      DAYRAN = DAYRAN+SUBRAN
      DAYRUN = DAYRUN+SUBRUN
      DAYSDH = DAYSDH+SUBSDH
      DAYSTP = DAYSTP+SUBSTP
      DAYTIM = DAYTIM+SUBTIM
      DAYTRA = DAYTRA+SUBTRA
      DAYUBC = DAYUBC+SUBUBC
      HOUR   = HOUR+DELSUB
      OLDRAT = RATE(NHR)
      DO I=1,NPT
        DAYSNK(I) = DAYSNK(I)+SUBSNK(I)
        DAYQL(I)  = DAYQL(I)+SUBQL(I)
        DAYQVH(I) = DAYQVH(I)+SUBQVH(I)
        DAYQVT(I) = DAYQVT(I)+SUBQVT(I)
      ENDDO
      IF (IHEAT .EQ. 1) THEN
        DAYHBE = DAYHBE+SUBHBE
        DAYLE  = DAYLE+SUBLE
        DAYRN  = DAYRN+SUBRN
        DAYSHF = DAYSHF+SUBSHF
        DAYSEN = DAYSEN+SUBSEN
        DAYQHW0 = DAYQHW0+SUBQHW0
        DO I=1,NPT
          DAYQHC(I) = DAYQHC(I)+SUBQHC(I)
          DAYQHV(I) = DAYQHV(I)+SUBQHV(I)
          DAYQHW(I) = DAYQHW(I)+SUBQHW(I)
        ENDDO
        IF (UPPERH .EQ. 0)                                              &
     &    CALL AIRTMP(TSTEPE,TAMEAN,TMAXLAST,TMINNEXT,TA)
      ENDIF
!
!     Write DELSUB summary to binary output file
!
      IF (NPRINT .EQ. 1) THEN
        IREC = IREC+1
        IF (IPLANT .EQ. 0) THEN
          IF (IHEAT .EQ. 0) THEN
            WRITE(LUB,REC=IREC) IDAY,HOUR,(H(J),THETA(J),SUBQL(J),      &
     &      SUBQVH(J),J=1,NPT),SUBINF,SUBRAN,SUBE,SUBTRA,SUBRUN,        &
     &      SPREV,TMOIST,SUBSTP,SUBPE,SUBPT,SUBTIM,SUBAST,SUBUBC,RHMEAN,&
     &      TAMEAN,HDRY,SUBRN,SUBSHF,SUBSEN,SHPREV,HSTORE,SUBQHW0,      &
     &      TA,TMAX,TMIN,VD_A,WIND,CLOUD,SR_MEAS,SUBLE,SUBSDH,SUBHBE
          ELSE
            WRITE(LUB,REC=IREC) IDAY,HOUR,(H(J),THETA(J),SUBQL(J),      &
     &      SUBQVH(J),SUBQVT(J),T(J),J=1,NPT),                          &
     &      SUBINF,SUBRAN,SUBE,SUBTRA,SUBRUN,                           &
     &      SPREV,TMOIST,SUBSTP,SUBPE,SUBPT,SUBTIM,SUBAST,SUBUBC,RHMEAN,&
     &      TAMEAN,HDRY,SUBRN,SUBSHF,SUBSEN,SHPREV,HSTORE,SUBQHW0,      &
     &      TA,TMAX,TMIN,VD_A,WIND,CLOUD,SR_MEAS,SUBLE,SUBSDH,SUBHBE,   &
     &      (SUBQHC(J),SUBQHW(J),SUBQHV(J),J=1,NPT)
          ENDIF
        ELSE
          WRITE(LUB,REC=IREC) IDAY,HOUR,(H(J),THETA(J),SUBQL(J),        &
     &    SUBQVH(J),SUBSNK(J),J=1,NPT),SUBINF,SUBRAN,SUBE,SUBTRA,SUBRUN,&
     &    SPREV,TMOIST,SUBSTP,SUBPE,SUBPT,SUBTIM,SUBAST,SUBUBC,RHMEAN,  &
     &    TAMEAN,HDRY,SUBSDH
        ENDIF
        IF (IDAY .EQ. IDEND .AND. N .EQ. NSTOP) THEN
          FINISH = .TRUE.
          GO TO 32673
        ENDIF
      ENDIF
!
!     Estimate new DELT to start the next DELSUB period
!
      DELT = MIN(1.0*DELMAX,(2.0*DELSUB)/SUBSTP)
      ENDDO
!----------------------------------------------------------------------
!     End of the NTOTAL Loop for this day
!
!     Add daily total to simulation total
!----------------------------------------------------------------------
32673 TAST   = TAST+DAYAST
      TSTP   = TSTP+DAYSTP
      TTIM   = TTIM+DAYTIM
      TPE    = TPE+DAYPE
      TE     = TE+DAYE
      IF (NGROW .GT. 0 .AND. NGROW .LE. NGDAYS) TETRAN = TETRAN+DAYE
      TIRR   = TIRR+DAYIRR
      TINF   = TINF+DAYINF
      TRUN   = TRUN+DAYRUN
      TTRA   = TTRA+DAYTRA
      TPT    = TPT+DAYPT
      TRAN   = TRAN+DAYRAN
      TPET   = TPET+PE+PT
      TSDH   = TSDH+DAYSDH
      TUBC   = TUBC+DAYUBC
      DO I=1,NPT
        TQL(I) = TQL(I)+DAYQL(I)
      ENDDO
      IF (IVAPOR .EQ. 1) THEN
        DO I=1,NPT
          TQVH(I) = TQVH(I)+DAYQVH(I)
          TQVT(I) = TQVT(I)+DAYQVT(I)
        ENDDO 
      ENDIF
      IF (IPLANT .EQ. 1) THEN
        DO I=1,NPT
          TSNK(I) = TSNK(I)+DAYSNK(I)
        ENDDO 
      ENDIF
      IF (IHEAT .EQ. 1) THEN
        THBEOLD = THBEOLD+DAYHBE
        TLE  = TLE+DAYLE
        TRN  = TRN+DAYRN
        TSHF = TSHF+DAYSHF
        TSEN = TSEN+DAYSEN
        TQHW0 = TQHW0+DAYQHW0
        DO I=1,NPT
          TQHC(I) = TQHC(I)+DAYQHC(I)
          TQHV(I) = TQHV(I)+DAYQHV(I)
          TQHW(I) = TQHW(I)+DAYQHW(I)
        ENDDO 
      ENDIF
!
!     Write daily summary to binary output file
!
      IREC = IREC+1
      IF (IPLANT .EQ. 0) THEN
        IF (IHEAT .EQ. 0) THEN
          WRITE(LUB,REC=IREC) IDAY,HOUR,(H(J),THETA(J),DAYQL(J),        &
     &    DAYQVH(J),J=1,NPT),DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,          &
     &    TPREV,TMOIST,DAYSTP,DAYPE,DAYPT,DAYTIM,DAYAST,DAYUBC,RHMEAN,  &
     &    TAMEAN,HDRY,DAYRN,DAYSHF,DAYSEN,DHPREV,HSTORE,DAYQHW0,        &
     &    TA,TMAX,TMIN,VD_A,WIND,CLOUD,SR_MEAS,DAYLE,DAYSDH,DAYHBE
        ELSE
          WRITE(LUB,REC=IREC) IDAY,HOUR,(H(J),THETA(J),DAYQL(J),        &
     &    DAYQVH(J),DAYQVT(J),T(J),J=1,NPT),                            &
     &    DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,                             &
     &    TPREV,TMOIST,DAYSTP,DAYPE,DAYPT,DAYTIM,DAYAST,DAYUBC,RHMEAN,  &
     &    TAMEAN,HDRY,DAYRN,DAYSHF,DAYSEN,DHPREV,HSTORE,DAYQHW0,        &
     &    TA,TMAX,TMIN,VD_A,WIND,CLOUD,SR_MEAS,DAYLE,DAYSDH,DAYHBE,     &
     &    (DAYQHC(J),DAYQHW(J),DAYQHV(J),J=1,NPT)
        ENDIF
      ELSE
          WRITE(LUB,REC=IREC) IDAY,HOUR,(H(J),THETA(J),DAYQL(J),        &
     &    DAYQVH(J),DAYSNK(J),J=1,NPT),DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,&
     &    TPREV,TMOIST,DAYSTP,DAYPE,DAYPT,DAYTIM,DAYAST,DAYUBC,RHMEAN,  &
     &    TAMEAN,HDRY,DAYSDH
      ENDIF
!
!     Write end-of-day status to screen 
!
      IF (ISCR .GT. 1) WRITE(*,9050) IYEAR,IDAY,DAYSTP
9050  FORMAT(' IYEAR =',I3,', IDAY ='I4,', Steps =',I5)
!
      IF (FINISH) GO TO 32706
      IF (IETOPT .EQ. 0) THEN
        IF (NWATER .GT. WINDEX .AND. IDAY .EQ. IRDAY) THEN
          CALL ZEROA(1,25,RTIME)
          CALL ZEROA(1,25,AMOUNT)
          READ(LUI) IRDAY,IRTYPE,EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
          WINDEX = WINDEX+1
        ENDIF
      ELSE
        IF (IDAY .EQ. IDEND) GO TO 32706
        CALL ZEROA(1,25,RTIME)
        CALL ZEROA(1,25,AMOUNT)
        TMAXLAST = TMAX
        READ(LUI) METDAY,TMAX,TMIN,TDEW,SR_MEAS,WIND,CLOUD,IRTYPE,      &
     &      EFICEN,NP,(RTIME(J),AMOUNT(J),J=1,NP)
        IF (METDAY .LT. IDEND) THEN
          READ(LUI) II,RR,TMINNEXT
          BACKSPACE 1
        ELSE
          TMINNEXT = TMIN
        ENDIF
        TAMEAN = (TMAX+TMIN)*.5
        PVD_A  = VD_A
        VD_A   = FUNVD(TDEW)
        SVD_A  = FUNVD(TAMEAN)
        RHMEAN = VD_A/SVD_A
        IF (ISHOPT .GT. 0 .AND. IHEAT .EQ. 0) THEN
          RRHA = RHA
          IF (ISHOPT .EQ. 1) RHA = RHMEAN
          IF (ISHOPT .EQ. 2) THEN
            TATIME = HOUR+0.5*DELSUB
            CALL AIRTMP(TATIME,TAMEAN,TMAXLAST,TMINNEXT,TA)
            RHA = VD_A/FUNVD(TA)
          ENDIF
          IF (RHA .GT. 1.0) RHA = 1.0
          IF (IEVOPT .EQ. 1) CALL HDRYCALC(RHA)
          IF (HDRY .LT. H(1)) THEN
            DELT = DELT*RAINIF
            IF (DELT .LT. DELMIN) DELT = DELMIN
          ENDIF
        ENDIF
        IF (IHEAT .EQ. 1 .AND. VD_A .GT. PVD_A) THEN
          DELT = DELT*RAINIF
          IF (DELT .LT. DELMIN) DELT = DELMIN
        ENDIF
        IF (NP .GT. 0) IRDAY = METDAY
      ENDIF
32706 ENDDO
!-----------------------------------------------------------------------
!     End of the Day Loop for this simulation
!
!     Write final simulation summary to binary output file
!-----------------------------------------------------------------------
      TERR = STORE+TINF-TE-TTRA-TQL(NPT)-TQVT(NPT)-TMOIST
      THBE = THPREV-HSTORE+TSHF+(TQHW0-TQHC(NPT)-TQHW(NPT)-TQHV(NPT))   &
     &       *10000.0
      IREC = IREC+1
      TTIM = TTIM/24.0
      WRITE(LUB,REC=IREC) IDAY,G,IPLANT,TPET,TPT,TTRA,TPE,TE,           &
     &    TETRAN,TRUN,TINF,TTIM,TRAN,APLIED,TIRR,TMOIST,TERR,           &
     &    TSTP,TAST,TUBC,TRN,TSHF,TSEN,TLE,THBE,TSDH,THPREV,HSTORE,     &
     &    (TQL(I),TQVH(I),TQVT(I),TQHC(I),TSNK(I),I=1,NPT),             &
     &    TQHW0,TQHW(NPT),HSOURCE,OSMPOT
!
!     If hysteresis is simulated, an ASCII output file of hysteresis
!     parameters is created, to be used as a restart file. The restart
!     file is given the same name as the input file but with the
!     extension "HRI", for Hysteresis Restart Information
!
      IF (IHYS .GT. 0) CALL HYSOUT
!
!     Write yearly simulation status to screen if requested
!
      IF (ISCR .GT. 0) WRITE(*,9100) IYS+IYEAR-1,TSTP
9100  FORMAT(' Completed year =',I5,', total steps =',I7)
      ENDDO
      WRITE(*,*) ' UNSAT-H simulation completed.'
      CLOSE(UNIT=LUB)
      CLOSE(UNIT=LUI,STATUS='DELETE')
999   STOP
!-----------------------------------------------------------------------
!     End of program UNSATH
!-----------------------------------------------------------------------
      END
      SUBROUTINE CLIP(ACHR)
!-----------------------------------------------------------------------
!     Processes filenames by clipping extraneous characters
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
!     End of subroutine CLIP supporting program UNSATH
!-----------------------------------------------------------------------
      END
      SUBROUTINE BLR(RH)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculates boundary layer resistance (RH [s/m]) to heat and vapor 
!     as a function of wind speed and measurement height, air 
!     temperature and measurement height, zero plane displacement, 
!     roughness lengths for head and momentum, and surface temperature
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (WIND .LE. 0.0) THEN
!     --when wind is zero, RH calculated as pure diffusion
!     --the factor of 10000 converts cm^2 to m^2
        RH = 10000.0*ZT/VAPDIF
        GO TO 20
      ENDIF
      IRHLIM = 5
      PSIM = 0.
      PSIH = 0.
      TSMEAN = 0.5*(TT(1)+T(1))
      DO I = 1,IRHLIM
        USTAR = VK*WIND/(RLNZU+PSIM)
        RH    = (RLNZT+PSIH)/(VK*USTAR)
        IF (I .EQ. 1) THEN
          RHO   = RH
          RHMIN = 0.1*RH
          RHMAX = 10.0*RH
        ELSE
          IF (RH .GT. RHMAX) THEN
            RH = RHMAX
            GO TO 20
          ELSE IF (RH .LT. RHMIN) THEN
            RH = RHMIN
            GO TO 20
          ENDIF
        ENDIF
        IF (I .LT. IRHLIM) THEN
          ZETA = -(VK*ZT*GRAV/(CHA*TA*USTAR*USTAR*USTAR))*              &
     &           (CHA*(TSMEAN-TA)/RH)
          IF (TSMEAN .LT. TA) THEN
            PSIH = 4.7*ZETA
            PSIM = PSIH
          ELSE
            PSIH = -2*LOG((1.0+SQRT(1.0-16.0*ZETA))/2.0)
            PSIM = 0.60*PSIH
          ENDIF
        ENDIF
      ENDDO
20    RETURN
!-----------------------------------------------------------------------
!     End of subroutine BLR
!-----------------------------------------------------------------------
      END
      SUBROUTINE DELCHK
!-----------------------------------------------------------------------
!     Called by UNSAT-H
!
!     Adjusts the time step according to the mass and heat balance error  
!     Specifically, it will increase or decrease the time step 
!     according to the ratio DMAXBA/DIFMAX
!
!     - if R = 0, the time step size is not changed
!     - if the time step size equals DELMIN and the solution
!       was unsuccessful (i.e., DIFMAX > DMAXBA), the program is
!       stopped rather than accepting the solution and moving on
!       (as was done in Version 2.0)
!     - if 0.9 < R < 1.0, R = 0.9; this change prevents the slow
!       decrease of time step reduction when DIFMAX is close to DMAXBA
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (DIFMAX .EQ. 0) THEN
        RWATER = 20.0
      ELSE
        RWATER = DMAXBA/ABS(DIFMAX)
      ENDIF
      R = RWATER
      IF (IHEAT .EQ. 1) THEN
        IF (DMAXHE .GT. 0 ) THEN
          IF (HEATBAL .EQ. 0) THEN
            RHEAT = 20.0
          ELSE
            RHEAT = DMAXHE/ABS(HEATBAL)
          ENDIF
          R = MIN(RWATER,RHEAT)
        ENDIF
      ENDIF
      IF (R .LT. 1.0) THEN
!-----------------------------------------------------------------------
!     Decrease the time step
!-----------------------------------------------------------------------
        IF (DELT .LE. DELMIN) THEN
          WRITE(*,*) ' Convergence criterion not met using DELMIN'
          WRITE(*,*) ' Simulation stopped in DELCHK'
          WRITE(*,6000) IDAY,ISWDIF,DMAXBA,DIFMAX,RWATER,DMAXHE,HEATBAL,&
     &      RHEAT,HOUR+TSUB-DELT,DELT,R
          STOP
        ENDIF
        RPRIM = R*R
        IF (RPRIM .GT. 0.9) RPRIM = 0.9
        IF (RPRIM .LT. 0.5) RPRIM = 0.5
        TSUB = TSUB - DELT
        DELT = DELT * RPRIM
        IF (DELT .LT. DELMIN) DELT = DELMIN
        IGO150 = 1
        TSUB = TSUB + DELT
        GO TO 1000
      ELSE
!-----------------------------------------------------------------------
!     Increase the time step
!-----------------------------------------------------------------------
        RPRIM = 0.5*(1+R)
        IF (RPRIM .GT. RFACT) RPRIM = RFACT
        DELT = DELT*RPRIM
        IF (DELT .GT. DELMAX) DELT = DELMAX
        TLEFT = DELSUB-TSUB
        IF (TLEFT .GT. 0. .AND. TLEFT .LE. 2.0*DELT+TIMRES) THEN
          IF (TLEFT .LE. DELT+TIMRES) THEN
            DELT = TLEFT
          ELSE
            DELT = TLEFT/2.0
          ENDIF
        ENDIF
      ENDIF
1000  RETURN
!-----------------------------------------------------------------------
!     Format statements
!-----------------------------------------------------------------------
6000  FORMAT(                                                           &
     &' IDAY   = ',I3,T24,   'ISWDIF  = ',I1,1P,/,                      &
     &' DMAXBA = ',E12.5,T24,'DIFMAX  = ',E12.5,T47,'RWATER = ',E12.5/, &
     &' DMAXHE = ',E12.5,T24,'HEATBAL = ',E12.5,T47,'RHEAT  = ',E12.5/, &
     &' TDONE  = ',E12.5,T24,'DELT    = ',E12.5,T47,'R      = ',E12.5) 
!-----------------------------------------------------------------------
!     End of subroutine DELCHK
!-----------------------------------------------------------------------
      END
      SUBROUTINE HDRYCALC(RHSURF)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calls RETENT, POLYKH
!
!     Calculates a value for HDRY based on the temperature of the soil
!     surface and the relative humidity of the atmosphere.
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      HDRY = -LOG(RHSURF)*TS/MGR
      H1   = H(1)
      H(1) = HDRY
      C1   = C(1)
      RKL1  = KL(1)
      THETA1 = THETA(1)
      CALL RETENT(1,1,H,1,THETA,1,C)
      CALL POLYKH(1,1,H,KL)
      DRYC  = C(1)
      DRYK  = KL(1)
      DRYTH = THETA(1)
      H(1)  = H1
      C(1)  = C1
      KL(1) = RKL1
      THETA(1) = THETA1
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine HDRYCALC
!-----------------------------------------------------------------------
      END
      SUBROUTINE NETRAD(IDAILY,INH)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculates net radiation (RN [W/m2]) 
!     as a function of latitude, time of year, time of day, 
!     air and surface temperature, atmospheric vapor density, albedo, 
!     cloud cover, surface and air emissivities, and extraterrestrial 
!     solar flux) in W/m2.
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (IDAILY .EQ. 1) THEN
        SINDEC = .3985*SIN(4.869+TPIY*METDAY+0.03345                    &
     &           *SIN(6.224+TPIY*METDAY))
        DEC = ASIN(SINDEC)
        COSDEC = COS(DEC)
        TANDEC = TAN(DEC)
        HS = ACOS(-TANLAT*TANDEC)
        SR_POT = DAYSEXT*(HS*SINLAT*SINDEC+COSLAT*COSDEC*SIN(HS))/PI
        SR_FRAC = SR_MEAS/SR_POT
!
!      Atmospheric emissivity calculated using equation from Campbell 
!      (1985) modified to accept vapor density units of g/cm3
!
        EMISS_A = 4.174*VD_A**(1.0/7.0)
        IF (ICLOUD .EQ. 0) THEN
          CLOUD = 2.330-3.330*SR_FRAC
          IF (CLOUD .GT. 1.0) CLOUD = 1.0
          IF (CLOUD .LT. 0.0) CLOUD = 0.0
        ENDIF
        IF (CLOUD .GT. 0) EMISS_A = (1.0-0.84*CLOUD)*EMISS_A+0.84*      &
     &      CLOUD
      ENDIF
      SR = SEXT*SR_FRAC*(SINLAT*SINDEC+COSLAT*COSDEC*                   &
     &     COS(PI*(TSTEPE-12.0)/12.0))
      IF (SR .LT. 0) SR = 0.0
!
!     ALBEDO and EMISS_S calculations are taken from van Bavel 
!     and Hillel (1976)
!
      AVETHE = (THETA(1)+TTHETA(1))/2.0
      IF (AVETHE .LT. 0.1) THEN
        ALBEDO = 0.25
      ELSE IF (AVETHE .GT. 0.25) THEN
        ALBEDO = 0.1
      ELSE
        ALBEDO = 0.1+0.25-AVETHE
      ENDIF
      EMISS_S = 0.9+0.18*AVETHE
!
!     These lines were used for the Tank project with Steve Simmons
! 
!     EMISS_A = EMISS_S+0.0351
!     EMISS_A = EMISS_S
      IF (EMISS_A .GT. 1.0) EMISS_A = 1.0
      IF (INH .EQ. 1) TAQ = TA*TA*TA*TA
      TSMEAN = 0.5*(TT(1)+T(1))
      TSQ = TSMEAN*TSMEAN*TSMEAN*TSMEAN
      RN = (1.0-ALBEDO)*SR+SB*(EMISS_A*TAQ-EMISS_S*TSQ)
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine NETRAD
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
      SUBROUTINE THERME(IBEGIN,IEND)
!-----------------------------------------------------------------------
!     Called by KVCALC
!
!     Calculates enhancement factor for soil vapor flow caused by
!     thermal gradients
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DO I=IBEGIN,IEND
        MATI = MAT(I)
        THETAI = THETA(I)
        EFA = EF(1,MATI)
        EFB = EF(2,MATI)
        EFC = EF(3,MATI)
        EFD = EF(4,MATI)
        EFE = EF(5,MATI)
        ENFACT(I) = EFA+EFB*THETAI-EFD*EXP(-1.0*(EFC*THETAI)**EFE)
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine THERME
!-----------------------------------------------------------------------
      END
      SUBROUTINE TRIDAG(IF,L,A,B,RC,RD,V)
!-----------------------------------------------------------------------
!     Called by HEATFLOW and UNSATH
!
!     Solves the tridiagonal solution matrix
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION A(*),B(*),RC(*),RD(*),V(*),BETA(M1),GAMM(M1)
      DATA BETA/M1*0/,GAMM/M1*0/
      BETA(IF)  = B(IF)
      GAMM(IF) = RD(IF)/BETA(IF)
      IFP1 = IF+1
      DO I=IFP1,L
        II = I - 1
        BETA(I) = B(I) - A(I)*RC(II)/BETA(II)
        GAMM(I) = (RD(I)-A(I)*GAMM(II))/BETA(I)
      ENDDO
      V(L) = GAMM(L)
      LAST = L - IF
      DO K = 1,LAST
        I = L - K
        V(I) = GAMM(I) - RC(I)*V(I+1)/BETA(I)
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine TRIDAG
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

      SUBROUTINE HEATFLOW(RH,INH)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     INH = 1  solves heat balance equations
!     INH = 2  calculates fluxes and heat balance error for current
!              time step solution
!     INH = 3  time step solution accepted; DELSUB values updated
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION A1(M1),A2(M1),A3(M1),ZY(M1)
      GO TO (100,200,300) INH
100   NPPT = NPT-1
      DEL = 1./DELSAV
!
!     Calculate coefficients for simulation of change in temperature
!
      DO II=2,NPPT
        I      = II-1
        III    = II+1
        Z31    = Z(III)-Z(I)
        KH1    = -KHMID(I)/(Z31*ZZ(I))
        KH2    = -KHMID(II)/(Z31*ZZ(II))
        AVCH   = 0.5*(CCH(II)+CH(II))
        A2(II) = AVCH*DEL-KH1-KH2
        A1(II) = KH1
        A3(II) = KH2
        TG1    = (TT(II)-TT(I))/ZZ(I)
        TG2    = (TT(III)-TT(II))/ZZ(II)
        ZY(II) = AVCH*TT(II)*DEL+HSOURCE*0.0036                         &
     &           -(-KKHMID(II)*TG2+KKHMID(I)*TG1)/Z31
        DELZ31 = 1/(DELSAV*Z31*2)
!
!     Add in the heat flux attributable to latent heat transfer
!
        IF (IVAPOR .EQ. 1) THEN
          DTHVDT = (THETAV(II)-TTHETAV(II))*DEL
          A1(II) = A1(II)-CHV*QV(I)*DELZ31
          A2(II) = A2(II)+CHV*0.5*DTHVDT+CHV*(QV(II)-QV(I))*DELZ31
          A3(II) = A3(II)+CHV*QV(II)*DELZ31
          ZY(II) = ZY(II)-(LHV0+CHV*(0.5*TT(II)-T0))*DTHVDT             &
     &             -((LHV0-CHV*T0)*(QV(II)-QV(I))                       &
     &             +(LHV0+CHV*(0.5*(TT(III)+TT(II))-T0))*QV(II)         &
     &             -(LHV0+CHV*(0.5*(TT(II)+TT(I))-T0))*QV(I))           &
     &             /(DELSAV*Z31)
        ENDIF
!
!     Add in the heat flux attributable to water flow
!
        IF (ICONVH .EQ. 1) THEN
          DTHDT  = (THETA(II)-TTHETA(II))*DEL
          A1(II) = A1(II)-CHW*QL(I)*DELZ31
          A2(II) = A2(II)+CHW*0.5*DTHDT                                 &
     &             +CHW*(QL(II)-QL(I))*DELZ31
          A3(II) = A3(II)+CHW*QL(II)*DELZ31
          ZY(II) = ZY(II)-CHW*(0.5*TT(II)-T0)*DTHDT                     &
     &             -(-CHW*T0*(QL(II)-QL(I))                             &
     &             +CHW*(0.5*(TT(III)+TT(II))-T0)*QL(II)                &
     &             -CHW*(0.5*(TT(II)+TT(I))-T0)*QL(I))/(DELSAV*Z31)
        ENDIF
      ENDDO
      IF (UPPERH .EQ. 0 .OR. UPPERH .EQ. 3) THEN
!
!     Surface node equation for heat flux b.c. (UPPERH = 0)
!     Surface fluxes such as net radiation, sensible heat, and latent
!     heat are calculated using the average time-step temperature,
!     Tavg = 0.5*(Tj + Tj-1). The surface flux units are W/m^2
!
        IF (UPPERH .EQ. 0) THEN
          SENFLX = CHA*(0.5*(TS+TT(1)-TA-TTA))/RH
          LE     = 0.0
!
!     When evaporation occurs, SSFLUX (cm) is converted into LE (W/m2)
!
          IF (SSFLUX .LT. 0) THEN 
            LE = -(LHV0+CHV*(0.5*(TS+TT(1))-T0))*SSFLUX/(DELSAV*0.36)
            SHFLUX = RN-SENFLX-LE
          ELSE
            SHFLUX = RN-SENFLX
          ENDIF
        ELSE
          SHFLUX = QHCTOP
        ENDIF
        Z1     = ZZ(1)
        KH2    = -KHMID(1)/(Z1*Z1)
        AVCH   = 0.5*(CH(1)+CCH(1))
        A1(1)  = 0.0
        A2(1)  = AVCH*DEL-KH2
        A3(1)  = KH2
        TG2    = (TT(2)-TT(1))/Z1
        ZY(1)  = AVCH*TT(1)*DEL+HSOURCE*0.0036                          &
     &           -(-KKHMID(1)*TG2-2.*SHFLUX*.36)/Z1
!
!     Add in the heat flux attributable to vapor flow
!
        IF (IVAPOR .EQ. 1) THEN 
          DTHVDT = (THETAV(1)-TTHETAV(1))*DEL
          A2(1)  = A2(1)+CHV*0.5*DTHVDT+CHV*QV(1)/(DELSAV*Z1*2)
          A3(1)  = A3(1)+CHV*QV(1)/(DELSAV*Z1*2)
          ZY(1)  = ZY(1)-(LHV0+CHV*(0.5*TT(1)-T0))*DTHVDT               &
     &             -((LHV0-CHV*T0)*QV(1)                                &
     &             +(LHV0+CHV*(0.5*(TT(2)+TT(1))-T0))*QV(1))            &
     &             /(DELSAV*Z1)
        ENDIF
!
!     Add in the heat flux attributable to water flow
!
        IF (ICONVH .EQ. 1) THEN
          DTHDT = (THETA(1)-TTHETA(1))*DEL
          A2(1) = A2(1)+CHW*0.5*DTHDT+CHW*QL(1)/(DELSAV*Z1*2)
          A3(1) = A3(1)+CHW*QL(1)/(DELSAV*Z1*2)
          ZY(1) = ZY(1)-CHW*(0.5*TT(1)-T0)*DTHDT-(CHW*QL(1)             &
     &            *(0.5*(TT(2)+TT(1))-2*T0))/(DELSAV*Z1)
          IF (SSFLUX .GT. 0) THEN
            RMOD = CHW*SSFLUX/(Z1*DELSAV)
            IF (UPPERH .EQ. 0) THEN
              ZY(1) = ZY(1)+2*RMOD*(0.5*(TA+TTA)-T0)
            ELSE IF (UPPERH .EQ. 3) THEN
              ZY(1) = ZY(1)+2*RMOD*(TSMEAN-T0)
            ELSE
              ZY(1) = ZY(1)+2*RMOD*(0.5*(TS+TT(1))-T0)
            ENDIF
          ENDIF
        ENDIF
        ISTART = 1
      ELSE 
!
!     Surface node equation for constant or specified temperature b.c.
!
        IF(UPPERH .EQ. 2) TS = TSMEAN+TSAMP*SIN(PI*(HOUR+TSUB-6.)/12.)
        ISTART = 2
        ZY(2)  = ZY(2)-A1(2)*TS
        A1(2)  = 0.
      ENDIF
!
!     LOWERH: 1) constant temperature gradient, 2) constant temperature,
!     and 3) constant flux
!
      IF (LOWERH .NE. 2) THEN
        ZN      = ZZ(NPPT)
        KH1     = -KHMID(NPPT)/(ZN*ZN)
        A3(NPT) = 0.0
        IF (LOWERH .EQ. 1) THEN
          QHLAST = -0.5*(KH(NPT)+KKHNPT)*TGRAD
        ELSE
          QHLAST = QHLEAK
        ENDIF
        AVCH    = 0.5*(CH(NPT)+CCH(NPT))
        A1(NPT) = KH1
        A2(NPT) = AVCH*DEL-KH1
        TG1     = (TT(NPT)-TT(NPPT))/ZN
        ZY(NPT) = AVCH*TT(NPT)*DEL+HSOURCE*0.0036                       &
     &            -(2*QHLAST+KKHMID(NPPT)*TG1)/ZN
!
!     Add in the heat flux attributable to vapor flow
!
        IF (IVAPOR .EQ. 1) THEN
          DTHVDT  = (THETAV(NPT)-TTHETAV(NPT))*DEL
          A1(NPT) = A1(NPT)-CHV*QV(NPPT)/(DELSAV*ZN*2)
          A2(NPT) = A2(NPT)+CHV*0.5*DTHVDT                              &
     &              +0.5*CHV*(QV(NPT)-QV(NPPT))/(DELSAV*ZN)
          ZY(NPT) = ZY(NPT)-(LHV0+CHV*(0.5*TT(NPT)-T0))*DTHVDT          &
     &             -((LHV0-CHV*T0)*(QV(NPT)-QV(NPPT))                   &
     &             +(LHV0+CHV*(0.5*(2*TT(NPT)+TGRAD*ZN)-T0))*QV(NPT)    &
     &             -(LHV0+CHV*(0.5*(TT(NPT)+TT(NPPT))-T0))*QV(NPPT))    &
     &             /(DELSAV*ZN)
        ENDIF
!
!     Add in the heat flux attributable to water flow
!
        IF (ICONVH .EQ. 1) THEN
          DTHDT   = (THETA(NPT)-TTHETA(NPT))*DEL
          A1(NPT) = A1(NPT)-CHW*QL(NPPT)/(DELSAV*ZN*2)
          A2(NPT) = A2(NPT)+CHW*0.5*DTHDT+CHW*(QL(NPT)                  &
     &              -0.5*QL(NPPT))/(DELSAV*ZN)
          ZY(NPT) = ZY(NPT)-CHW*(0.5*TT(NPT)-T0)*DTHDT                  &
     &              -(-CHW*T0*(QL(NPT)-QL(NPPT))                        &
     &              +CHW*(TT(NPT)-T0)*QL(NPT)                           &
     &              -CHW*(0.5*(TT(NPT)+TT(NPPT))-T0)*QL(NPPT))          &
     &              /(DELSAV*ZN)
        ENDIF
        NPPT = NPT
      ELSE
        ZY(NPPT) = ZY(NPPT)-A3(NPPT)*T(NPT)
        A3(NPPT) = 0.0
      ENDIF
!
!     Solution of tridiagonal matrix
!
      CALL TRIDAG(ISTART,NPPT,A1,A2,A3,ZY,T)
!-----------------------------------------------------------------------
!     Solution estimate complete for this iteration
!-----------------------------------------------------------------------
!     T limits are -40 C (-40 F) lower and 90 C (194 F) upper.
!     If T goes outside these limits, the time step is reduced.
!     Future option: add convergence criteria for temperature?
!
      DO I=1,NPT
        IF (T(I) .LT. 233 .OR. T(I) .GT. 363) THEN
          IF (DELT .GT. DELMIN) THEN
            TSUB = TSUB-DELT
            DELT = MAX(0.5*DELT,DELMIN)
            TSUB = TSUB+DELT
            IGO150 = 1
          ELSE
            WRITE(*,*) ' T limits exceeded and DELT < DELMIN'
            WRITE(*,*) ' Simulation stopped in HEATFLOW' 
            STOP
          ENDIF
        ENDIF
      ENDDO
!-----------------------------------------------------------------------
!     End of the iteration loop 
!-----------------------------------------------------------------------
      GO TO 900
!-----------------------------------------------------------------------
!     INH = 2
!     Calculate fluxes and heat balance error for this time step
!       heat flux across the surface (SHFLUX [W/m2])
!       heat fluxes within the soil (QHcvw [J/cm2/h])
!-----------------------------------------------------------------------
200   DO II=2,NPT
        I   = II-1
        ZI  = ZZ(I)
        TG  = (T(II)-T(I))/ZI
        TTG = (TT(II)-TT(I))/ZI
        QHC(I) = 0.5*(-KHMID(I)*TG-KKHMID(I)*TTG)
      ENDDO
!
!     Calculate the heat flux attributable to vapor flow
!
      IF (IVAPOR .EQ. 1) THEN
        DO I=1,NPT-1
          QHV(I) = (LHV0+CHV*(0.25*(T(I+1)+T(I)+TT(I+1)+TT(I))-T0))     &
     &             *QV(I)/DELSAV
        ENDDO
        QHV(NPT) = (LHV0+CHV*(0.5*(T(NPT)+TT(NPT)+TGRAD*ZZ(NPT-1))-T0)) &
     &             *QV(NPT)/DELSAV
      ENDIF
!
!     Calculate the heat flux attributable to liquid water flow
!
      IF (ICONVH .EQ. 1) THEN
        QHW0 = 0
        IF (SSFLUX .GT. 0) THEN
          IF (UPPERH .EQ. 0) THEN
            AVGTS = 0.5*(TA+TTA)
          ELSE IF (UPPERH .EQ. 3) THEN
            AVGTS = TSMEAN
          ELSE
            AVGTS = 0.5*(TS+TT(1))
          ENDIF
          QHW0 = CHW*(AVGTS-T0)*SSFLUX/DELSAV
        ENDIF
        DO I=1,NPT-1
          QHW(I) = CHW*(0.25*(T(I+1)+T(I)+TT(I+1)+TT(I))-T0)            &
     &             *QL(I)/DELSAV
        ENDDO
        QHW(NPT) = CHW*(0.5*(T(NPT)+TT(NPT))-T0)*QL(NPT)/DELSAV
      ENDIF
!
!     Calculate surface heat flux (SHFLUX) for specified surface
!     temperature conditions. SHFLUX does not include convective flux.
!
      IF (UPPERH .EQ. 1 .OR. UPPERH .EQ. 2) THEN
        SHFLUX = HSOURCE*THICK(1)*0.01+(QHC(1)                          &
     &          +(0.5*(CH(1)+CCH(1))*(TS-TT(1))*THICK(1))/DELSAV)/0.36
        IF (IVAPOR .EQ. 1)                                              &
     &    SHFLUX = SHFLUX+(QHV(1)+((LHV0+CHV*(0.5*(TS+TT(1))-T0))       &
     &             *(THETAV(1)-TTHETAV(1))*THICK(1))/DELSAV)/0.36
        IF (ICONVH .EQ. 1)                                              &
     &    SHFLUX = SHFLUX+(QHW(1)-QHW0+(CHW*(0.5*(TS+TT(1))-T0)         &
     &             *(THETA(1)-TTHETA(1))*THICK(1))/DELSAV)/0.36
      ENDIF
!
!     Calculate bottom conductive heat flux, QHC(NPT), [J/cm2/h]
!
      IF (LOWERH .EQ. 1) THEN
        QHC(NPT) = 0.5*(-KH(NPT)-KKHNPT)*TGRAD
      ELSE IF (LOWERH .EQ. 2) THEN
        QHC(NPT) = QHC(NPT-1)+HSOURCE*THICK(NPT)*0.0036
        THCKTM = THICK(NPT)/DELSAV
        IF (IVAPOR .EQ. 1) QHC(NPT) = QHC(NPT)+QHV(NPT-1)-QHV(NPT)      &
     &    -(LHV0+CHV*(0.5*(T(NPT)+TT(NPT))-T0))                         &
     &    *(THETAV(NPT)-TTHETAV(NPT))*THCKTM
        IF (ICONVH .EQ. 1) QHC(NPT) = QHC(NPT)-QHW(NPT)+QHW(NPT-1)      &
     &    -CHW*(0.5*(T(NPT)+TT(NPT))-T0)*(THETA(NPT)-TTHETA(NPT))       &
     &    *THCKTM
      ELSE
        QHC(NPT) = QHLEAK
      ENDIF
!
!     Estimate heat storage in the profile (SHSTORE) and the heat
!     balance error (HEATBAL) in J/m2.
!
      SHSTORE = 0.0
      HGEN    = 0.0
      DO I=1,NPT
        SHSTORE = SHSTORE+CH(I)*(T(I)-T0)*THICK(I)
        HGEN = HGEN+HSOURCE*THICK(I)
        IF (IVAPOR .EQ. 1) SHSTORE = SHSTORE+LHV0*THETAV(I)*THICK(I)
      ENDDO
      HGEN = HGEN*36.0
      SHSTORE = SHSTORE*10000.0
      HEATBAL = (SHFLUX*3600.0+10000.0*(HGEN+QHW0-QHW(NPT)-QHV(NPT)     &
     &          -QHC(NPT)))*DELSAV+HSTORE-SHSTORE
      DHSTORE = 0.0
      DO I=1,NPT
        DHSTORE = DHSTORE+0.5*(CH(I)+CCH(I))*(T(I)-TT(I))*THICK(I)
        IF (IVAPOR .EQ. 1)                                              &
     &    DHSTORE = DHSTORE+(LHV0+CHV*(0.5*(T(I)+TT(I))-T0))            &
     &    *(THETAV(I)-TTHETAV(I))*THICK(I)
        IF (ICONVH .EQ. 1)                                              &
     &    DHSTORE = DHSTORE+CHW*(0.5*(T(I)+TT(I))-T0)                   &
     &    *(THETA(I)-TTHETA(I))*THICK(I)
      ENDDO
      HEATBL = (SHFLUX*0.36+HGEN+QHW0                                   &
     &          -QHW(NPT)-QHV(NPT)-QHC(NPT))*DELSAV-DHSTORE
!-----------------------------------------------------------------------
!     End of convergence loop for HEATFLOW.
!-----------------------------------------------------------------------
      GO TO 900
!-----------------------------------------------------------------------
!     INH = 3
!     Time step solution is accepted, DELSUB totals are updated
!     The surface heat fluxes are converted from [W/m2] to [J/m2],
!       other fluxes are [J/cm2]
!     Old heat flow parameters are reset to new values
!-----------------------------------------------------------------------
300   SUBHBE  = SUBHBE+HEATBAL
      SUBLE   = SUBLE+LE*3600*DELSAV
      SUBRN   = SUBRN+RN*3600*DELSAV
      SUBSHF  = SUBSHF+SHFLUX*3600*DELSAV
      SUBSEN  = SUBSEN+SENFLX*3600*DELSAV
      SUBQHW0 = SUBQHW0+QHW0*DELSAV
      DO I=1,NPT
        SUBQHC(I) = SUBQHC(I)+QHC(I)*DELSAV
        SUBQHV(I) = SUBQHV(I)+QHV(I)*DELSAV
        SUBQHW(I) = SUBQHW(I)+QHW(I)*DELSAV
      ENDDO
900   RETURN
!-----------------------------------------------------------------------
!     End of subroutine HEATFLOW
!-----------------------------------------------------------------------
      END
      SUBROUTINE AIRTMP(TATIME,RTAMEAN,TMAXLAST,TMINNEXT,RTA)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculates the air temperature as a function of time of day.
!     Assumes air temperature varies sinusoidally during each day with
!     the minimum air temperature occurring at 0300 and maximum air
!     temperature occurring at 1500 h
!     TATIME   = time of day (h) at which to evaluate air temperature
!     RTA      = air temperature (C)
!     TMAX     = maximum air temperature of the current day
!     TMAXLAST = maximum air temperature of the previous day
!     TMIN     = minimum air temperature of the current day
!     TMINNEXT = minimum air temperature of the next day
!     RTAMEAN  = mean air temperature of current day (TMAX+TMIN)*.5
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      IF (TATIME .LT. 3.) THEN
        TMEANE = (TMAXLAST+TMIN)*0.5
        IF (TMIN .GT. TMAXLAST) THEN
          RTA = TMEANE+(TMIN-TMEANE)*COS((3.0-TATIME)*PI/12.0)
        ELSE
          RTA = TMEANE+(TMAXLAST-TMEANE)*COS((15.0-TATIME)*PI/12.0)
        ENDIF
      ELSE IF (TATIME .GE. 15.) THEN
        TMEANL = (TMAX+TMINNEXT)*0.5
        IF (TMINNEXT .GT. TMAX) THEN
          RTA = TMEANL+(TMINNEXT-TMEANL)*COS((3.0-TATIME)*PI/12.0)
        ELSE
          RTA = TMEANL+(TMAX-TMEANL)*COS((15.0-TATIME)*PI/12.0)
        ENDIF
      ELSE
        RTA = RTAMEAN+(TMAX-RTAMEAN)*COS((15.0-TATIME)*PI/12.0)
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine AIRTMP
!-----------------------------------------------------------------------
      END
      SUBROUTINE ZEROA(IBEGIN,IEND,ARRAY)
!-----------------------------------------------------------------------
!     Called by UNSAT-H
!
!     Sets ARRAY to zero
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION ARRAY(M1)
      DO I=IBEGIN,IEND
        ARRAY(I) = 0.0
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine ZEROA
!-----------------------------------------------------------------------
      END
      SUBROUTINE ZEROI(K1,K2,K3,K4)
!-----------------------------------------------------------------------
!     Called by UNSAT-H
!
!     Sets INTEGER parameters to zero
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      K1 = 0.0
      K2 = 0.0
      K3 = 0.0
      K4 = 0.0
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine ZEROI
!-----------------------------------------------------------------------
      END
      SUBROUTINE ZEROR(R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13,R14)
!-----------------------------------------------------------------------
!     Called by UNSAT-H
!
!     Sets REAL parameters to zero
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      R1  = 0.0
      R2  = 0.0
      R3  = 0.0
      R4  = 0.0
      R5  = 0.0
      R6  = 0.0
      R7  = 0.0
      R8  = 0.0
      R9  = 0.0
      R10 = 0.0
      R11 = 0.0
      R12 = 0.0
      R13 = 0.0
      R14 = 0.0
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine ZEROR
!-----------------------------------------------------------------------
      END
      SUBROUTINE ROOTS(NGROW,MDEPTH)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculates maximum rooting depth and the root-density function
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      I = 1
400   IF (NGROW .GT. NTROOT(I) .AND. I .LT. MXROOT) THEN
        I = I+1
        GO TO 400
      ENDIF
      MDEPTH = I
      IF (MDEPTH .GT. 1) THEN
!
!     Estimate root density factor (RDF) at each node of root depth
!     Note that roots are not allowed to grow in node 1 (surface node)
!
        RTOT = 0.
        DO I=2,MDEPTH
          RDF(I) = RLD(Z(I))
          RTOT   = RTOT+RDF(I)*THICK(I)
        ENDDO
        DO I=1,MDEPTH
          RDF(I) = RDF(I)/RTOT
        ENDDO
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine ROOTS
!-----------------------------------------------------------------------
      END
      SUBROUTINE KVCALC(IBEGIN,IEND,TSOIL,VC)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculate vapor conductivities
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      FUNVD(TD) = EXP(46.440973-(6790.4985/TD)-6.02808*LOG(TD))
      IF (IVAPOR .EQ. 1) THEN
        IF (IHEAT .EQ. 0) THEN
!
!     VC was modified for isothermal vapor flow (i.e., already
!     multiplied by SVD*MGR/TSOIL because TSOIL is constant) 
!
          DO I=IBEGIN,IEND
            AIRCON = DMAX1(0.0D0,THET(MAT(I))-THETA(I)-ENTAIR(I))
            CALL RELHUM(H(I),TSOIL,RHUMID)
            KV(I) = VC*AIRCON*RHUMID
          ENDDO
        ELSE
          CALL THERME(IBEGIN,IEND)
          DO I=IBEGIN,IEND
            TI = T(I)
            AIRCON = DMAX1(0.0D0,THET(MAT(I))-THETA(I)-ENTAIR(I))
            CALL RELHUM(H(I),TI,RHUMID)
            SVD    = FUNVD(TI)
            DSVDDT = SVD*(-6.02808+6790.4985/TI)/TI
            KV(I)  = VC*SVD*AIRCON*RHUMID*MGR/TI
            KVT(I) = VC*AIRCON*ENFACT(I)*RHUMID*DSVDDT
          ENDDO
        ENDIF
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine KVCALC
!-----------------------------------------------------------------------
      END
      SUBROUTINE INTERK(IBEGIN,IEND,KEST,UP,DOWN)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     1) Computes internodal conductivity, KLMID, based on KEST options
!     2) Computes internodal conductivities, KVMID and KVTMID, based
!        on KVEST options (currently, KVEST set equal to KEST)
!     3) Computes internodal thermal conductivity, KHMID, using
!        arithmetic averaging
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      LOGICAL LMDIFF
      MM = 0
      DO I=IBEGIN,IEND
!
!     Compute KLMID 
!
      RKL1 = KL(I)
      RKL2 = KL(I+1)
      X1   = H(I)
      X2   = H(I+1)
      IF (KEST .GT. 3) THEN
        M = MAT(I)
        IF (M .EQ. MM) THEN
          LMDIFF = .FALSE.
        ELSE
          LMDIFF = .TRUE.
          MM = M
        ENDIF
      ENDIF
!
!     Arithmetic mean of KL, KEST=1
!
      IF (KEST .EQ. 1) THEN
        KVEST = 1
        IF (UP .EQ. 0.5) THEN
          KLMID(I) = (RKL1+RKL2)*0.5
        ELSE IF (-(X1+Z(I)) .GT. -(X2+Z(I+1))) THEN
          KLMID(I) = RKL1*UP+RKL2*DOWN
        ELSE
          KLMID(I) = RKL1*DOWN+RKL2*UP
      ENDIF
!
!     Harmonic mean of KL, KEST=2
!
      ELSE IF (KEST .EQ. 2) THEN
        KVEST    = 2
        KLMID(I) = 2*(RKL1*RKL2/(RKL1+RKL2))
!
!     Geometric mean of KL, KEST=3
!
      ELSE IF (KEST .EQ. 3) THEN
        KVEST    = 3
        KLMID(I) = (RKL1*RKL2)**0.5
        IF (KLMID(I) .EQ. 0) KLMID(I) = EXP((LOG(RKL1)+LOG(RKL2))/2)
      ENDIF
!
!     Compute KVMID and KVTMID 
!
      IF (IVAPOR .EQ. 1) THEN
        RKV1 = KV(I)
        RKV2 = KV(I+1)
        IF (IHEAT .EQ. 1) THEN
          RKVT1 = KVT(I)
          RKVT2 = KVT(I+1)
        ENDIF
!
!     Arithmetic mean, KVEST=1
!
        IF (KVEST .EQ. 1) THEN
          IF (UP .EQ. 0.5) THEN
            KVMID(I) = (RKV1+RKV2)*0.5
            IF (IHEAT .EQ. 1) KVTMID(I) = (RKVT1+RKVT2)*0.5
          ELSE IF (-(X1+Z(I)) .GT. -(X2+Z(I+1))) THEN
            KVMID(I) = RKV1*UP+RKV2*DOWN
            IF (IHEAT .EQ. 1) KVTMID(I) = RKVT1*UP+RKVT2*DOWN
          ELSE 
            KVMID(I) = RKV1*DOWN+RKV2*UP
            IF (IHEAT .EQ. 1) KVTMID(I) = RKVT1*DOWN+RKVT2*UP
          ENDIF
!
!     Harmonic mean, KVEST=2
!
        ELSE IF (KVEST .EQ. 2) THEN
            IF (RKV1 .EQ. 0.0 .AND. RKV2 .EQ. 0.0) THEN
              KVMID(I) = 0.0
            ELSE
              KVMID(I) = 2*(RKV1*RKV2/(RKV1+RKV2))
            ENDIF
            IF (IHEAT .EQ. 1) THEN
              IF (RKVT1 .EQ. 0.0 .AND. RKVT2 .EQ. 0.0) THEN
                KVTMID(I) = 0.0
              ELSE
                KVTMID(I) = 2*(RKVT1*RKVT2/(RKVT1+RKVT2))
              ENDIF
            ENDIF
!
!     Geometric mean, KVEST=3
!
          ELSE IF (KVEST .EQ. 3) THEN
            KVMID(I) = (RKV1*RKV2)**0.5
            IF (IHEAT .EQ. 1) KVTMID(I) = (RKVT1*RKVT2)**0.5
          ENDIF
        ENDIF 
!
!     Compute internodal thermal conductivities [J/cm/hr/K]
!
        IF (IHEAT .EQ. 1) KHMID(I) = (KH(I)+KH(I+1))*0.5
      ENDDO
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine INTERK
!-----------------------------------------------------------------------
      END
      SUBROUTINE RESET(IBEGIN,IEND,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,  &
     &    A12,A13,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11,R12,R13)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Copy R values into A.  Used to
!       a) initialize parameters at start of iteration loop
!       b) reset parameters after a time step reduction
!
!     Parameter Number (array A or R)
!      1  C
!      2  H
!      3  THETA
!      4  KLMID
!      5  KVMID
!      6  KVTMID
!      7  T
!      8  CH
!      9  KHMID
!     10  THETAV
!     11  KL(NPT)
!     12  KH(NPT)
!     13  KVT(NPT)
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION A1(M1),A2(M1),A3(M1),A4(M1),A5(M1),A6(M1),A7(M1),A8(M1),&
     &    A9(M1),A10(M1),R1(M1),R2(M1),R3(M1),R4(M1),R5(M1),R6(M1),     &
     &    R7(M1),R8(M1),R9(M1),R10(M1)
      DO I=IBEGIN,IEND
        A1(I) = R1(I)
        A2(I) = R2(I)
        A3(I) = R3(I)
        A4(I) = R4(I)
      ENDDO
      IF (IVAPOR .EQ. 1) THEN
        DO I=IBEGIN,IEND
          A5(I) = R5(I)
        ENDDO
      ENDIF
      IF (IHEAT .EQ. 1) THEN
        DO I=IBEGIN,IEND
          A6(I) = R6(I)
          A7(I) = R7(I)
          A8(I) = R8(I)
          A9(I) = R9(I)
        ENDDO
        IF (IVAPOR .EQ. 1) THEN
          DO I=IBEGIN,IEND
            A10(I) = R10(I)
          ENDDO
        ENDIF
      ENDIF
      A11 = R11
      A12 = R12
      A13 = R13
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine RESET
!-----------------------------------------------------------------------
      END
      SUBROUTINE TRSINK(MDEPTH,PTSUB,NSINK,SINK,TRANSP)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculates sink term at each node; totals sink terms to get
!     transpiration
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      DIMENSION SINK(*)
      TRANSP = 0.0
      CALL ZEROA(1,MDEPTH,SINK)
!
!     Transpiration is assumed to be zero during irrigation or rainfall
!
      IF (PTSUB .EQ. 0.0) NSINK = 0
!
!     Estimate sink term by a modified FEDDES method
!
      IF (NSINK .GT. 0) THEN
        DO I=1,MDEPTH
          MATI = MAT(I)
          TH   = (THETA(I)+TTHETA(I))*0.5
          THW  = THETAW(MATI)
          THD  = THETAD(MATI)
          THN  = THETAN(MATI)
          IF (TH .GT. THW .AND. TH .LE. THN) THEN
            IF (TH .LT. THD) THEN
              ALPHAF = (TH-THW)/(THD-THW)
            ELSE
              ALPHAF = 1.0
            ENDIF
          ELSE
            ALPHAF = 0.0
          ENDIF
          SINK(I) = PTSUB*ALPHAF*RDF(I)
          TRANSP  = TRANSP+SINK(I)*THICK(I)
        ENDDO
      ENDIF
      RETURN
!-----------------------------------------------------------------------
!     End of subroutine TRSINK
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
      SUBROUTINE FLUX(RH)
!-----------------------------------------------------------------------
!     Called by UNSATH
!
!     Calculates liquid and vapor water fluxes between nodes
!-----------------------------------------------------------------------
      INCLUDE 'unsath.inc'
      CALL ZEROA(1,NPT,QL)
      CALL ZEROA(1,NPT,QT)
      IF (IVAPOR .EQ. 1) THEN
        CALL ZEROA(1,NPT,QV)
        CALL ZEROA(1,NPT,QVH)
        IF (IHEAT .EQ. 1) CALL ZEROA(1,NPT,QVT)
      ENDIF
      DO II=2,NPT
        I     = II-1
        ZI    = ZZ(I)
        PG    = (H(II)-H(I))/ZI
        PPG   = (HH(II)-HH(I))/ZI
        QL(I) = DELSAV*0.5*(KLMID(I)*(PG+G)+KKLMID(I)*(PPG+G))
        IF (IVAPOR .EQ. 1) THEN
          QVH(I) = DELSAV*0.5*(KVMID(I)*PG+KKVMID(I)*PPG)
          IF (IHEAT .EQ. 1) THEN
          TG  = (T(II)-T(I))/ZI
          TTG = (TT(II)-TT(I))/ZI
          QVT(I) = -DELSAV*0.5*(KVTMID(I)*TG+KKVTMID(I)*TTG)
          ENDIF
          QV(I) = QVH(I)+QVT(I)
        ENDIF
        QT(I) = QL(I)+QV(I)
      ENDDO
!
!     SSFLUX is surface infiltration (cm) for the time step.
!     Estimate from THETA and geometry if saturated.
!     If the infiltration was from rainfall, however, do not let it
!     exceed the actual rainfall (i.e. ISAVUP .NE. 1)
!
      IF (UPPER .EQ. 1) THEN
        SSFLUX = QT(1)+(THETA(1)-TTHETA(1))*THICK(1)
        IF (ISAVUP .NE. 1) THEN
          IF (H(1) .LE. HIRRI) THEN
            SSFLUX = MIN(SSFLUX,SFLUX*DELSAV)
          ELSE
            IF (IHEAT .EQ. 0 .AND. ISATUK .EQ. 0) THEN
              SSFLUX = MAX(SSFLUX,SFLUX*DELSAV)
              SSFLUX = MIN(0.0D0,SSFLUX)
            ENDIF
          ENDIF
        ENDIF
      ELSE
        IF (ISATUK .EQ. 0) THEN
          IF (IEVOPT .EQ. 3) THEN
            RHS = EXP(-H(1)*MGR/TS)
            IF (RHA  .LT. 1.0) EFLUX = -PESUB*(RHS-RHA)/(1-RHA)
            IF (EFLUX .GT. 0.0) EFLUX = 0.0
            SFLUX = 0.5*(EEFLUX+EFLUX)
          ENDIF
          IF (IHEAT .EQ. 1 .AND. UPPERH .EQ. 0) THEN
            SFLUX = 3.6D5*(VD_A-((VVD_S+VD_S)/2.0))/(WATDEN*RH)
            IF (SFLUX .GT. 0.0) SFLUX = 0.0
          ENDIF
        ENDIF
        SSFLUX = SFLUX*DELSAV
      ENDIF
!
!     Deep drainage flux
!
      IF (LOWER .EQ. 2 .OR. LOWER .EQ. 5) THEN
        QL(NPT) = QT(NPT-1)
      ELSE
        QL(NPT) = QLAST*DELSAV
      ENDIF
      QT(NPT) = QL(NPT)
      IF (IHEAT .EQ. 1 .AND. LOWERH .EQ. 1) THEN
        IF (IVAPOR .EQ. 1) THEN
          QVT(NPT) = -DELSAV*0.5*(KVT(NPT)+KKVTNPT)*TGRAD
          QV(NPT)  = QVT(NPT)
          QT(NPT)  = QT(NPT)+QV(NPT)
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
!     End of subroutine FLUX
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
