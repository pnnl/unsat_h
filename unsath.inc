!-----------------------------------------------------------------------
!     File UNSATH.INC
!       Dimensions
!       Logical unit parameters
!       Record length specification
!       Common blocks
!       Equivalences
!     Program Units:  all
!-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
!-----------------------------------------------------------------------
!     PARAMETER assignments
!-----------------------------------------------------------------------
!
!     M1 is the maximum dimension for nodal arrays
!
      PARAMETER (M1 = 250)
!
!     M2 is the maximum number of soil types
!
      PARAMETER (M2 = 10)
!
!     M3 is the maximum number of polynomials for any soil type.
!     A minimum value of 1 is needed for all non-polynomial options.
!
      PARAMETER (M3 = 4)
!
!     M4 is the maximum number of coefficients per soil function
!       Must be > or = 5 for Haverkamp, BC, and VG functions
!       Must be > or = 7 for Rossi-Nimmo functions, mod. BC function
!       Must be > or = 14 for modified VG function
!       Must be > or = MAXPOL*MAXCOE for the polynomial option
!
      PARAMETER (M4 = 20)
!
!     IRLMOD modifies the output file record length.
!       IRLMOD = 1  VAX and PC single precision
!       IRLMOD = 2  VAX and PC double precision
!       IRLMOD = 8  SUN Workstation
!
      PARAMETER (IRLMOD = 2)
!
!     Logical unit assignments
!       LUI = binary input file unit (input to UNSATH)
!       LUB = binary output file unit from DATAINH;  
!             binary output file unit from UNSATH
!       LUR = read (interactive input unit)
!       LUS = screen (interactive and additional output unit)
!       LUW = output files from DATAOUT
!       LUH = output file containing hysteresis restart information
!       LUX = scratch file
!
      PARAMETER (LUI = 1, LUB = 2, LUR = 5, LUS = 6, LUW = 7, LUH = 8)
      PARAMETER (LUX = 3)
!
!     Internally calculated parameters (used in equivalence statements)
!       M2M3 = maximum allowed number of polynomials for all soil types
!       M2M4 = maximum number of soil coefficients per function for
!              all soil types
!       M5 = total number of water-related soil coefficients for SWPA
!       M6 = array dimension for TCON and EF
!       M7 = total number of heat-related soil coefficients for SHPA
!       M8 = total number of values in two hysteresis arrays
!       M9 = total number of PET/MET and PRECIP input file names
!
      PARAMETER (M2M3 = M2*M3, M2M4 = M2*M4, M5 = 9*M2+4*M2M3+2*M2M4)
      PARAMETER (M6 = 5*M2, M7 = 2*M1+11*M2, M8 = 7*M1, M9 = 100)
!-----------------------------------------------------------------------
!     End of PARAMETER assignments
!
!     Begin data type assignments
!-----------------------------------------------------------------------
      CHARACTER*80 IFILE,OFILE,BIFILE,RFILE,TMFILE,HYFILE,PMFN(M9),     &
     &  PRFN(M9),TITLE,DUMMY,STIME*11,SDATE*11
      DIMENSION AIRINT(M2),THET(M2),THTR(M2),XDIVTH(M2M3),CREGTH(M2M4), &
     &  AIRINK(M2),SK(M2),XDIVKH(M2M3),CREGKH(M2M4),THTA(M1),HM(M1),    &
     &  TCON(5,M2),CHS(M2),EF(5,M2),CHSOIL(M1),ENFACT(M1),AMOUNT(25),   &
     &  RTIME(25),NSUBTH(M2),NDEGTH(M2M3),NSUBKH(M2),NDEGKH(M2M3)
      REAL*8 LOG10E,MGR,LAT,LHV0,LHV,LE,KKLMID(M1),KKLNPT,KL(M1),       &
     &  KKHNPT,KVT(M1),KKVTNPT,KVTMID(M1),KH(M1),KKVTMID(M1),KH1,KH2,   &
     &  KHMID(M1),KKHMID(M1),KV(M1),KVMID(M1),KKVMID(M1),KLMID(M1)
      INTEGER DAYSTP,TSTP,DAYAST,TAST,DAYUBC,TUBC,                      &
     &  UPPER,UPPERH,DAYSDH,TSDH
      LOGICAL LUSLUW
!-----------------------------------------------------------------------
!     End data type assignments
!
!     Begin COMMON and EQUIVALENCE assignments
!-----------------------------------------------------------------------
!
!     Character variables
!
      COMMON /CHAR/ IFILE,OFILE,BIFILE,RFILE,TMFILE,HYFILE,PMFN,PRFN,   &
     &  TITLE,DUMMY,STIME,SDATE
!
!     Meteorological and surface variables
!
      COMMON /CMET/RMDATA(8,366),PET(366),FPET(24),ALT,PMB,ALBEDO,VK,   &
     &  PI,CHW,ZH,ZM,ZT,ZU,D,CHA,SB,SEXT,DAYSEXT,TAQ,GRAV,HPR,          &
     &  SINDEC,COSDEC,LAT,SINLAT,COSLAT,TANLAT,TPIY,LHV,                &
     &  RN,SUBSEN,SUBSHF,SUBRN,SUBLE,SUBHBE,VD_A,SR_FRAC,               &
     &  TMAX,TMIN,TMEAN,TA,TTA,SR_MEAS,WIND,CLOUD,EMISS_A,WATDEN,       &
     &  TSMEAN,TSAMP,RHMEAN,QHCTOP,QHLEAK,SENFLX,                       &
     &  LHV0,T0,CHV,LE,RLNZU,RLNZT,ICLOUD,NWATER,IRAIN,IETOPT,UPPER,    &
     &  UPPERH,ISHOPT
      EQUIVALENCE (T(1),TS)
!
!     Soil property variables; soil water property array (SWPA)
!
      COMMON /CSOIL/SWPA(M5),HDRY,HIRRI,SATK,DRYK,SATC,DRYC,SATTH,      &
     &  DRYTH,VAPDIF,WCOMPR,RKMOD,KOPT,IVAPOR,MAXPOL,MAXCOE,MAT(M1),INC
      EQUIVALENCE (SWPA(1),NSUBTH),                                     &
     &  (SWPA(1+  M2),AIRINT),                                          &
     &  (SWPA(1+2*M2),THET),                                            &
     &  (SWPA(1+3*M2),THTR),                                            &
     &  (SWPA(1+4*M2),THTA),                                            &
     &  (SWPA(1+5*M2),HM),                                              &
     &  (SWPA(1+6*M2),NDEGTH),                                          &
     &  (SWPA(1+6*M2+  M2M3),XDIVTH),                                   &
     &  (SWPA(1+6*M2+2*M2M3),CREGTH),                                   &
     &  (SWPA(1+6*M2+2*M2M3+M2M4),NSUBKH),                              &
     &  (SWPA(1+7*M2+2*M2M3+M2M4),AIRINK),                              &
     &  (SWPA(1+8*M2+2*M2M3+M2M4),SK),                                  &
     &  (SWPA(1+9*M2+2*M2M3+M2M4),NDEGKH),                              &
     &  (SWPA(1+9*M2+3*M2M3+M2M4),XDIVKH),                              &
     &  (SWPA(1+9*M2+4*M2M3+M2M4),CREGKH)
!
!     Parameters for modified van Genuchten function (Fayer-Simmons 95)
!
      COMMON /CMODVG/VGA,VGN,VGM,P,CK1,CK2,CK3,HMM,WHM,HO,WHO,XO,HXO,   &
     &  HC,WHC,GAMMA,GAMMAM,GAMMAX,RLAMB,SR,WO,GSXO,RI3BWO,SA
!
!     Time parameters
!
      COMMON /CTIM/AMOUNT,RTIME,DELT,DIFMAX,DELSAV,DELSUB,TSUB,TSTART,  &
     &  DMAXBA,DMAXHE,DELMAX,DELMIN,RFACT,HTIME,HOUR,OUTTIM,HEATBAL,    &
     &  DAYTIM,STOPHR,TIMRES,TSTEPB,TSTEPM,TSTEPE,NDAYS,NDAY,DAYSTP,    &
     &  TSTP,IDAY,DAYAST,TAST,DAYSDH,TSDH,DAYUBC,TUBC,IDBEG,IDEND,      &
     &  METDAY,N,ISWDIF,IGO150,IN,NTOTAL,NYEARS,IYS
!
!     Water variables
!
      COMMON /CWATER/H(M1),NTROOT(M1),HH(M1),RDF(M1),THETA(M1),         &
     &  TTHETA(M1),C(M1),CC(M1),T(M1),TT(M1),QT(M1),PEVAPO(366),        &
     &  THICK(M1),ZZ(M1),QL(M1),QVH(M1),QVT(M1),QV(M1),THETAV(M1),      &
     &  TTHETAV(M1),QWLEAK(366),                                        &
     &  DAYQL(M1),DAYQVH(M1),DAYQVT(M1),TQL(M1),TQVH(M1),               &
     &  TQVT(M1),BBHEAD(50,2),SFLUX,SSFLUX,TPREV,TMOIST
!
!     Heat variables
!
      COMMON /CHEAT/SHPA(M7),KVT,KVTMID,KH,KKVTMID,KHMID,KKHMID,        &
     &  KL,KV,KVMID,KKVMID,KLMID,KKLMID,CH(M1),CCH(M1),                 &
     &  QHC(M1),SUBQHC(M1),DAYQHC(M1),TQHC(M1),                         &
     &  QHV(M1),SUBQHV(M1),DAYQHV(M1),TQHV(M1),                         &
     &  QHW(M1),SUBQHW(M1),DAYQHW(M1),TQHW(M1),KKVTNPT,                 &
     &  KKHNPT,KH1,KH2,KKLNPT,DHPREV,HSTORE,THPREV,SHSTORE,IHEAT,ICONVH
      EQUIVALENCE (SHPA(1),CHSOIL),                                     &
     &  (SHPA(1+  M1),ENFACT),                                          &
     &  (SHPA(1+2*M1),CHS),                                             &
     &  (SHPA(1+2*M1+  M2),TCON),                                       &
     &  (SHPA(1+2*M1+6*M2),EF)
!
!     Plant variables
!
      COMMON /CPLANT/ PTRANS(366),DLAI(366),DAYSNK(M1),TSNK(M1),        &
     &  THETAW(M2),THETAD(M2),THETAN(M2),PETPC(5),AA,B1,B2,DAYTRA,      &
     &  BIOMAS,BARE,IPLANT,MXROOT,NFPET,NSOW,NHRVST,NGDAYS
!
!     Hysteresis variables
!
      COMMON/CHYSTR/RASW(M1,7),RHSW(M1,7),RSWAW(M1),SARW(M1),ENTAIR(M1),&
     &  ASWMXH(M1),SARWA(M2),RAWA(M2),ALFACT(M2),HYSHPH(M2),            &
     &  ALPHAD,ALPHAI,XN,XM,SM,PHI,XXM,RSW,SRAW,HYSMXH,CWW,PERMW,       &
     &  HYSTOL,AIRTOL,OMTOL3,OPTOL3,ASW,ESW,SATW,SW,ASAT,RAW,SARWI      &
      COMMON/CHYSTI/JJH(M1),JH(M1),IPSW(M1),IPATHA(M2),                 &
     &  IHYS,LSW,IPSW1,JJJH,IPATH
!
!     Real variables
!
      COMMON /CREAL/LOG10E,MGR,QHLAST,TGRAD,QLAST,PESUB,                &
     &  EEFLUX,EFLUX,RHA,VVD_S,VD_S,SVD_S,HSOURCE,OSMPOT,               &
     &  DAYINF,DAYRAN,DAYE,DAYRUN,PE,PT,DOV,DAYRN,DAYSHF,DAYSEN,        &
     &  DAYHBE,DAYLE,SUBQHW0,DAYQHW0,TQHW0
!
!     Geometry variables
!
      COMMON /CGEOM/Z(M1),G,NPT
!
!     Integer variables
!
      COMMON /CINTEG/LOWER,ISAVUP,ISATUK,IEVOPT,LOWERH,ILEAP,INDEXB
!
!     Output processing variables
!
      COMMON /COUT/NODE(7),LOUT,NPRINT,MAXREC,IREC,NDEPTH,NCHR,IERROR
!-----------------------------------------------------------------------
!     End of file UNSATH.INC
!-----------------------------------------------------------------------
