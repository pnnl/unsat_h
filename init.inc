!-----------------------------------------------------------------------
!     File init.inc
!     Program Units: ALL
!     All common blocks an constants are initialized, with the
!     exception of some constants set in DATAINH.F
!-----------------------------------------------------------------------
      DATA DELT,DIFMAX,DELSAV,ISWDIF,DELSUB,TSUB,TSTART,                &
     &  IDAY,N,IGO150,RFACT,HTIME,HOUR,OUTTIM /14*0/
      DATA RTIME/25*0/,RRHA,RRHS/2*0.5/,RH/10.0/,AMOUNT/25*0/
      DATA Z/M1*0/,H/M1*0/,NTROOT/M1*0/,HH/M1*0/,RDF/M1*0/,             &
     &  THETA/M1*0/,TTHETA/M1*0/,C/M1*0/,CC/M1*0/,T/M1*0/,              &
     &  TT/M1*0/, KL/M1*0/
      DATA FPET/24*0/,THETAW/M2*0/,THETAD/M2*0/,THETAN/M2*0/,           &
     &  KKHMID/M1*0/,PTRANS/366*0/,PEVAPO/366*0/,QWLEAK/366*0/
      DATA QT/M1*0/,CCH/M1*0/,QV/M1*0/,THICK/M1*0/,ZZ/M1*0/,            &
     &  QHC/M1*0/,SUBQHC/M1*0/,DAYQHC/M1*0/,TQHC/M1*0/,QL/M1*0/,        &
     &  QHV/M1*0/,SUBQHV/M1*0/,DAYQHV/M1*0/,TQHV/M1*0/,                 &
     &  QHW/M1*0/,SUBQHW/M1*0/,DAYQHW/M1*0/,TQHW/M1*0/,                 &
     &  TTHETAV/M1*0/,THETAV/M1*0/
      DATA RASW/M8*0/,RHSW/M8*0/,RSWAW/M1*0/,SARW/M1*0/,ENTAIR/M1*0/,   &
     &  ASWMXH/M1*0/,SARWA/M2*0/,RAWA/M2*0/,                            &
     &  ALFACT/M2*0/,HYSHPH/M2*0/
      DATA DMAXBA,DMAXHE,DELMAX,DELMIN,HIRRI,HDRY,TMOIST,LOG10E,MGR,    &
     &  SATK,DRYK,SATC,DRYC,SATTH,DRYTH,AA,B1,B2,SSFLUX,QHLAST,TGRAD,   &
     &  RKMOD,VAPDIF,HSTORI,HSTORE,HSOURCE,OSMPOT /27*0.0/
      DATA IPLANT,LOWER,NDAYS,NDAY,KOPT,IVAPOR,MATN,NPT,                &
     &  NGRAV,MAXPOL,MAXCOE,IEVOPT,INC,MXROOT,NWATER,                   &
     &  NPRINT,IETOPT,ISHOPT,IHEAT,LOWERH,UPPERH,ICLOUD,ICONVH,         &
     &  IDBEG,IDEND /25*0/
      DATA START,FINISH /2*.TRUE./
      DATA MAT /M1*0/,TCON/M6*0/,CHS/M2*0/,                             &
     &  EF/M6*0/,CHSOIL/M1*0/,ENFACT/M1*0/,CH/M1*0/,KH/M1*0/
      DATA SWPA /M5*0/
      DATA SHPA /M7*0/
      DATA DAYMET,ALT,LAT,PMB,ALBEDO,ZH,ZM,ZT,ZU,D,DAYSEXT,TAQ,         &
     &  SINDEC,COSDEC,SINLAT,COSLAT,TANLAT,TPIY,LHV,                    &
     &  RN,SUBSEN,SUBSHF,SUBRN,SUBLE,SUBHBE,VD_A,SR_FRAC,               &
     &  TMAX,TMIN,TMEAN,TA,SR_MEAS,WIND,CLOUD,EMISS_A,                  &
     &  TSMEAN,TSAMP,KKHNPT,ILEAP,RHMEAN,METDAY,QHCTOP,QHLEAK,          &
     &  SENFLX,LE /45* 0.0/
      DATA PET /366*0/
      DATA DAYSTP,TSTP,LOUT,QHW0,SUBQHW0,DAYQHW0,TQHW0,                 &
     &  NTOTAL,MAXREC,IREC,IERROR,NCHR,DAYAST,TAST,                     &
     &  DAYUBC,TUBC,UPPER,NDEPTH,DAYSDH,TSDH /20*0/
      DATA DAYINF,DAYRAN,DAYE,DAYTRA,DAYRUN,TPREV,                      &
     &  PE,PT,DAYTIM,STOPHR,DAYRN,DAYSHF,DAYSEN,                        &
     &  DAYHBE,DAYLE/15*0.0/
      DATA PI/3.14159265358979/
!     The following parameters are fixed within UNSAT-H:
!     CHA = Volumetric specific heat of air     [J/m3/K]
!     CHV = Volumetric specific heat of vapor   [J/cm3/K]
!     CHW = Volumetric specific heat of water   [J/cm3/K]
!     GASCON = Gas constant                     [erg/mole/K]
!     GRAV = Gravity constant                   [m/s2]
!     MOLAR = Molar mass of water               [g/mole]
!     LHV0 = Latent heat of vaporization at     [J/g]
!            reference temperature T0
!     SB = Stephan-Boltzman constant            [W/m2/K4]
!     SEXT = Solar constant                     [W/m2]
!     T0 = Reference temperature                [K]
!     VK = von Karmens constant                 [-]
!     WATDEN = Density of water                 [g/cm3]
      DATA CHA    /1200.0/
      DATA CHV    /1.8/
      DATA CHW    /4.18/
      DATA GASCON /8.3143E7/
      DATA GRAV   /9.807/
      DATA MOLAR  /18.0/
      DATA LHV0   /2500./
      DATA SB     /5.67E-8/
      DATA SEXT   /1360.0/
      DATA T0     /273.16/
      DATA VK     /0.4/
      DATA WATDEN /1.0/
!
!     WCOMPR is the compressibility of water, 1/cm
!
      DATA WCOMPR /4.58D-08/
!
!     Hourly precip. rate when using only daily precip. data, cm/h
!      
      DATA HPR /1.0/
!
!     Hysteresis tolerances (min. air content, cm3/cm3;
!     head change for reversal, cm)
!
      DATA AIRTOL /0.00001/, HYSTOL/0.1/
!
!     Minimum time resolution (hour)
!
      DATA TIMRES /1.0D-12/
!
!     Output logical - if user writes to file, warning messages sent to
!     screen and file
!
      LUSLUW = .FALSE.
!-----------------------------------------------------------------------
!     End of file init.inc
!-----------------------------------------------------------------------
