DEFINE_CALL GLA(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                OUTPUT)
DEFINE_CALL PENMAN(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,   OUTPUT,OUTPUT)
DEFINE_CALL EVAPTR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,INPUT,        OUTPUT,OUTPUT,OUTPUT)
DEFINE_CALL DRUNIR(INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT,INPUT, ...
                   INPUT,INPUT,INPUT,              OUTPUT,OUTPUT,OUTPUT)

TITLE LINTUL2 CASSAVA VERSION 1.0 by Ezui K.S. & P.A. Leffelaar, 2016
*----------------------------------------------------------------------*
*  LINTUL, Light INTerception and UtiLization simulator                *
*          A simple general crop growth model, which simulates dry     *
*          matter production as the result of light interception and   *
*          utilization with a constant light use efficiency.           *
*  LINTUL2 is an extended version of LINTUL1 (the version of LINTUL    *
*          for optimal growing conditions). LINTUL2 includes a simple  *
*          water balance for studying effects of drought. The water    *
*          balance can be found in section 6 of the program, and the   *
*          effect of drought on light use efficiency in section 4.     *
*                                                                      *
*          Example for spring wheat                                    *
*                                                                      *
*  DLO-Research Institute for Agrobiology and Soil Fertility (AB-DLO)  *
*  Dept of Theor. Prod. Ecology, Wageningen Agric. Univ. (TPE-WAU)     *
*                                                                      *
*  Reference: Spitters, C.J.T. & A.H.C.M. Schapendonk, 1990.           *
*  Evaluation of breeding strategies for drought tolerance in potato   *
*  by means of crop growth simulation. Plant and Soil 123: 193-203.    *
*                                                                      *
*  Cite this model as follows:                                         *
*                                                                      *
*  Ezui, K.S., P.A. Leffelaar, A.C. Franke, A. Mando & K.E. Giller,    *
*  2017. Understanding water-limited yield of cassava in southern Togo.*
*  In: Ezui, K.S. (ed), Understanding the productivity of cassava in   *
*  West Africa. PhD Thesis, Wageningen University, NL, pp. 108-148.    *
*----------------------------------------------------------------------*

***   1. Initial conditions and run control

INITIAL

*     Initial conditions
INCON ZERO = 0.; ROOTDI = 0.1
      WCUTTINGIP = WCUTTINGUNIT * NCUTTINGS
      WLVI = WCUTTINGIP * FLV_CUTT
      LAII = WLVI * SLAI
      WAI  = 1000. * ROOTDI * WCI

SET EMERG = 0.

* WCUTTINGIP is the stem cutting weight at planting
* NCUTTINGS is the number of cuttings planted per m2
* WCUTTINGUNIT is the average weight per cutting (g)
* NCUTTINGS and WCUTTING help calculate initial stem weight (WSTI)

*     Run control
FINISH TSUM > 4320.
TIMER STTIME = 113.; FINTIM = 435.; DELT = 1.; PRDEL = 1.

* FINISH TSUM is determined assuming a 12-month life cycle and a daily TSUM of 12 deg days (opt temp 27oC - TBASE 15oC)
* Thus, FINISH TSUM = 12 deg days x 360 days. The 360 days are obtained as 12 months x 30 days/month

TRANSLATION_GENERAL DRIVER='EUDRIV'
PRINT DAP, TSUM, TSUMCROP, LAI, WSO, WGTOTAL, TRANRF, DORMANCY, PTRAN, TRAIN, DAVTMP

*PRINT DAP, TSUM, TSUMCROP, EMERG, LAI, GLV, WLV, WLVG, WLVD, WST, WSO, WRT, WGTOTAL, CUMUL_PARINT, TRANRF, DORMANCY, ...
*      PUSHDORMREC, PUSHREDIST, PUSHREDISTEND, WSOREDISTFRAC, PUSHREDISTSUM, WLVGREDISTCUMUL, RDRDV, RDRSH, RDRSD, ...
*      RROOTD, ROOTD, REDISTMAINTLOSSCUMUL, FRACWCUTTING, FRACWSO, FRACWST, FRACWLV, FRACWRT, FRACWTOT, ...
*      WC, WCWP, WCCR, WCSD, PTRAN, RAIN, TRAIN, DAVTMP, WCUTTING


***   2. Environmental data and temperature sum

DYNAMIC
** Total cum evapotranspiration reference values
CUMEVP = INTGRL(ZERO,TOTEVP)
TOTEVP = PEVAP + PTRAN

** Total cum evapotranspiration actual values
CUMEVA = INTGRL(ZERO,TOTEVA)
TOTEVA=  EVAP + TRAN
WP     = WSO/NOTNUL(CUMEVA)
CUMUL_EVAP = INTGRL(ZERO, EVAP)
CUMUL_TRAN = INTGRL(ZERO, TRAN)

WEATHER WTRDIR='C:\LINTUL\LINTUL1\Weather\';CNTR='SVK';ISTN=6;IYEAR=2013

*     Reading weather data from weather file:
*     RDD    Daily global radiation        J/(m2*d)
*     TMMN   Daily minimum temperature     degree C
*     TMMX   Daily maximum temperature     degree C
*     VP     Vapour pressure               kPa
*     WN     Wind speed                    m/s
*     RAIN   Precipitation                 mm

      DTR    = RDD/1.E+6
      TRAIN  = INTGRL(ZERO,RAIN)
      TDTR   = INTGRL(ZERO,DTR)
      DAVTMP = 0.5 * (TMMN + TMMX)
      DTEFF  = MAX ( 0., DAVTMP-TBASE )
      TSUM   = INTGRL(ZERO, RTSUM)
      RTSUM  = DTEFF*INSW(TIME-DOYPL, 0., 1.)
      TSUMCROP   = INTGRL(ZERO, RTSUMCROP)
      RTSUMCROP = DTEFF * EMERG

* The accumulation of TSUM for the whole growth starts at planting to induce TSUM accumulation for emergence (EMERGTSUM).
* But another TSUM accumulation is used from emergence to harvest (TSUMCROP)

* The model assumes EMERG (Emergence) occurs when TSUM accumulation for EMERG (EMERGTSUM) has reached the
* optimal emergence TSUM requirement (OPTEMERGTSUM) and WC>WCWP, then stays ON for the rest of the growth process

      EMERGTSUM = INTGRL(ZERO, REMERGTSUM)

      REMERGTSUM = DTEFF * INSW(EMERG - 1.,1.,0.)*INSW(TIME - DOYPL,0.,1.)

      DAP = TIME - DOYPL

      EMERGFAVOURABLE = REAAND(TSUM - OPTEMERGTSUM, WC-WCWP) * INSW(-WST, 0.,1.)

* EMERGFAVOURABLE is defined as intermediary push function to indicate the appropriate conditions in which EMERG should occur.
* EMERGFAVOURABLE is used to define the ZEROCONDITION in the EVENT function below
* At emergence (EMERG = 1.), initial weight of the stem cutting, the stem, the roots, the leaves and storage roots are defined.
* The stem cutting is used as the main source of DM for the stem, the roots and the leaves. There is no storage roots production
* at emergence (FSO_CUTT = 0. at EMERG)

EVENT
   ZEROCONDITION 1. - EMERGFAVOURABLE
   NEWVALUE  EMERG = 1.
   NEWVALUE  WCUTTING = WCUTTING *(1. - FST_CUTT - FRT_CUTT - FLV_CUTT - FSO_CUTT)
   NEWVALUE  WST = WST + WCUTTINGIP * FST_CUTT
   NEWVALUE  WRT = WRT + WCUTTINGIP * FRT_CUTT
   NEWVALUE  WLVG = WLVG + WCUTTINGIP * FLV_CUTT
   NEWVALUE  WSO = WSO + WCUTTINGIP * FSO_CUTT
ENDEVENT


***   3. Stem cutting growth (weight decline)

      WCUTTING = INTGRL(WCUTTINGIP, RWCUTTING)
      RWCUTTING = -RDRWCUTTING * WCUTTING * INSW(WCUTTING - WCUTTINGMIN, 0.,1.) * TRANRF * EMERG * (1. - DORMANCY)
      WCUTTINGMIN = WCUTTINGMINPRO * WCUTTINGIP

* Stem cutting weight declines from its initial value at emergence to a min weight (WCUTTINGMIN) from which
* no more weight decline occurs. This min weight corresponds in proportion to WCUTTINGMINPRO, a min fraction of WCUTTING
* from which there is no more cutting weight decline.
* The amount of DM declined from WCUTTINGIP is added to DM produced from photosynthesis (GTOTAL)

***   4. Leaf growth and senescence

** Specific leaf area

      SLA = SLA_MAX * FRACSLACROPAGE
      FRACSLACROPAGE = AFGEN(FRACSLATB, TSUMCROP)

* According to Keating et al., 1982, and as confirmed during our field experiment, SLA was low in the early stage of the crop's
* growth but increased as the crop age increased: 57% of the max SLA at 4MAP, 65% at 8MAP and 100% at 11MAP.
* These fractions were used to define FRACSLACROP in relation to the crop physiological time from emergence.
* SLA_MAX is the maximum SLA achieved in later stage of the crop's development.

** Leaf area growth and senescence

      CALL GLA(DTEFF,TSUMCROP,LAII,RGRL,DELT,SLA,LAI,...
               GLV,TSUMLA_MIN,TRANRF,WC,WCWP,...
               RWCUTTING,FLV,LAIEXPOEND,DORMANCY,...
               GLAI)


      GLV   = FLV * (GTOTAL + ABS(RWCUTTING)) + RREDISTLVG * PUSHREDIST

      DLV   = (WLVG * RDR - RWSOFASTRANSLSO) * (1. - DORMANCY)

      RWSOFASTRANSLSO = WLVG * RDR * FASTRANSLSO * (1. - DORMANCY)

* The model assumes that before shedding, leaves (to be dropped) reallocate part of their content to storage roots, following
* Johnson et al. 1986, as implemented in SUBSTOR-Potato.
* RWSOFASTRANSLSO is the rate of DM reallocation from leaves to storage roots before leaf dropping from the plant at senescence.
* FASTRANSLSO is the fraction of the senesced leaf DM that is reallocated to storage roots.
* The model assumes that FASTRANSLSO applies to all leaf death conditions (age, shade and enhanced shedding due to severe drought).

      DLAI  = LAI * RDR * (1. - FASTRANSLSO) * (1. - DORMANCY)


** Adding the effect of leaf shedding due to drought (RDRSD)
* Leaf death is assumed to be caused either by leaf age, shade or severe drought leading to enhanced shedding.

      RDR   = MAX(RDRDV, RDRSH, RDRSD) * INSW(TSUMCROPLEAFAGE - TSUMLLIFE, 0.,1.)
      TSUMCROPLEAFAGE = INTGRL(ZERO, RTSUMCROPLEAFAGE)
      RTSUMCROPLEAFAGE = DTEFF * EMERG - (TSUMCROPLEAFAGE/DELT) * PUSHREDIST

      RDRSD = RDRB * ENHSHED
      RDRDV = INSW(TSUMCROPLEAFAGE-TSUMLLIFE, 0., AFGEN(RDRT, DAVTMP))
      RDRSH = LIMIT(0., RDRSHM, RDRSHM * (LAI-LAICR) / LAICR)

* RDRB is a parameter defined as basic relative death rate of leaves

      ENHSHED = MAX(INSW(WC-WCSD, 1.,0.), INSW(WC-WCWET,0.,1.))*INSW(TSUMCROPLEAFAGE-FRACTLLFENHSH*TSUMLLIFE, 0.,1.)

* WCSD is WC at severe drought.
* The model assumes that leaf shedding enhancement (ENHSHED) will occur after a certain period of time because leaf life
*(from emergence to abscission) is between 40 and 210 days (Alves, 2002). In our case, we used TSUMLLIFE for this assumed period.
* TSUMLLIFE is the developmental time from emergence of a new leaf to its shedding.
* A separate TSUM (TSUMCROPLEAFAGE) is used to check leaf age beyond which leaf senescence should start.
* The model assumes that leaf shedding can be enhanced from half life to full life duration
* using FRACTLLFENHSH (The fraction of the leaf life from which enhanced shedding can start. Its values: 0.5 - 1).
* TSUMSBR is the developmental time for the occurrence of the first branch (sympodial branch) that marks the start of
* the reproductive stage (Matthews and Hunt, 1994).

      RLAI  = GLAI - DLAI
      LAI   = INTGRL(ZERO, RLAI)

***   5. Light interception and total crop growth rate

      PARINT = 0.5 * DTR * (1. - EXP(-K_ext*LAI))

      GTOTAL = LUE * PARINT * TRANRF * (1. - DORMANCY)

      LUE = LUE_OPT * AFGEN(TTB, DAVTMP)

      CUMUL_PARINT = INTGRL(ZERO, PARINT)

* Cassava growth is optimal within a 25 and 29 deg C (Alves, 2002). Within that temperature range, the model assumes
* an optimal LUE value (LUE_OPT). We assume that the maximum LUE values achieved during our field experiment (> 75th percentile)
* were reached within that optimal temperature range.

***   6. Growth rates and dry matter production of plant organs

      FRTWET = AFGEN( FRTTB, TSUMCROP )
         FRTMOD = MAX( 1., 1./(TRANRF+0.5) )
      FRT    = FRTWET * FRTMOD
         FSHMOD = (1.-FRT) / (1.-FRT/FRTMOD)
      FLV    = AFGEN( FLVTB, TSUMCROP ) * FSHMOD
      FST    = AFGEN( FSTTB, TSUMCROP ) * FSHMOD
      FSO    = AFGEN( FSOTB, TSUMCROP ) * FSHMOD
      FTOT   = FRT + FLV + FST + FSO

      RWST   = (GTOTAL + ABS(RWCUTTING)) * FST
      RWLVG  = (GTOTAL + ABS(RWCUTTING)) * FLV - DLV + RREDISTLVG * PUSHREDIST
      RWRT   = (GTOTAL + ABS(RWCUTTING)) * FRT
      RWSO   = (GTOTAL + ABS(RWCUTTING)) * FSO + RWSOFASTRANSLSO - RREDISTSO

* The rate of dry matter reallocation from leaves before shedding (RWSOFASTRANSLSO) is added to storage roots DM growth rate (RWSO)
* The rate of dry matter redistribution for leaf production at the recovery from dormancy (RREDISTSO) is substrated from RWSO.
* The rate of production of new leaves at the recovery from dormancy (RREDISTLVG) is added to RWLVG, provided PUSHREDIST is active.

      WLVG   = INTGRL( ZERO, RWLVG)
      WLVD   = INTGRL( ZERO, DLV  )
      WST    = INTGRL( ZERO, RWST )
      WSO    = INTGRL( ZERO, RWSO )
      WSOFASTRANSLSO = INTGRL( ZERO, RWSOFASTRANSLSO )
      WRT    = INTGRL( ZERO, RWRT )
      WLV = WLVG + WLVD
      WTOPS = WLV + WST + WCUTTING
      WTOPSGREEN = WLVG + WST + WCUTTING
      WGTOTALHARV = WTOPSGREEN + WSO

      WGTOTAL = WTOPS + WSO + WRT

* REDISTMAINTLOSSCUMUL is the amount of dry matter lost during the redistribution from storage roots to leave

         WSOTHA = WSO / 100.
         WSTTHA = WST / 100.
         WLVGTHA = WLVG / 100.
         WLVDTHA = WLVD / 100.
         WLVTHA = WLV / 100.
         WRTTHA = WRT / 100.
         WTOPSTHA = WTOPS / 100.
         WTOPSGREENTHA = WTOPSGREEN / 100.
         WGTOTALHARVTHA = WGTOTALHARV / 100.
         WGTOTALTHA = WGTOTAL / 100.

* Harvest index (HI)
      HI = WSO / NOTNUL(WGTOTALHARV)

* Checking of biomass proportions in each plant organ
      FRACWCUTTING = WCUTTING / NOTNUL(WGTOTAL)
      FRACWSO = WSO / NOTNUL(WGTOTAL)
      FRACWST = WST / NOTNUL(WGTOTAL)
      FRACWRT = WRT / NOTNUL(WGTOTAL)
      FRACWLV = WLV / NOTNUL(WGTOTAL)
      FRACWLVG = WLVG / NOTNUL(WGTOTAL)
      FRACWLVD = WLVD / NOTNUL(WGTOTAL)
      FRACWTOT = FRACWCUTTING + FRACWSO + FRACWST + FRACWRT + FRACWLV

** Root depth

      RROOTD = RRDMAX * INSW( WC-WCWP, 0., 1. ) * INSW(ROOTD-ROOTDM, 1.,0.) * EMERG

      ROOTD  = INTGRL( ROOTDI, RROOTD)

***   7. Soil moisture balance

      EXPLOR = 1000. * RROOTD * WCFC
      RNINTC = MIN( RAIN, FRACRNINTC*LAI )

      CALL PENMAN( DAVTMP,VP,DTR,LAI,WN,RNINTC, ...
                   PEVAP,PTRAN)
      CALL EVAPTR( PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,WCFC,WCWET,WCST,...
                   TRANCO,DELT,WCSD,...
                   WCCR,EVAP,TRAN)
          TRANRF = TRAN / NOTNUL(PTRAN)
      CALL DRUNIR( RAIN,RNINTC,EVAP,TRAN,IRRIGF,...
                   DRATE,DELT,WA,ROOTD,WCFC,WCST,...
                   DRAIN,RUNOFF,IRRIG)
      RWA = (RAIN+EXPLOR+IRRIG) - (RNINTC+RUNOFF+TRAN+EVAP+DRAIN)
      WA  = INTGRL( WAI,RWA)
      WC  = 0.001 * WA/ROOTD
      WCSD = WCWP * TWCSD

* TWCSD is a parameter accounting for how much higher WCSD is, relatively to WCWP

**** 7b. Dormancy and recovery from dormancy

      DORMANCY = MAX(REANOR(WC-WCSD, LAI-LAI_MIN), PUSHDORMREC)*(1. - PUSHREDIST)*INSW(TSUMCROP - TSUMSBR, 0.,1.)

* DORMTSUM is TSUM accumulation during the dormancy
* Just after the dormancy, DORMTSUM is cleared or re-initialised to zero with the activation of PUSHREDIST

      RDORMTSUM = DTEFF * DORMANCY - (DORMTSUM/DELT) * PUSHREDIST

      DORMTSUM = INTGRL(ZERO, RDORMTSUM)

* Duration of the dormancy
      DURDORM = INTGRL(ZERO, RDORMTIME)
      RDORMTIME = DORMANCY

* The recovery from dormancy is preceded by WC back to WCCR, which activates PUSHDORMREC and
* TSUM accumulation until the delay period (DELREDIST) is reached. When DELREDIST is reached, this triggers DM distribution
* from storage roots to leaves

      PUSHDORMREC = REAAND(WC-RECOV*WCCR,WC-WCWP)*INSW(-DORMTSUM, 1.,0.)*(1. - PUSHREDIST)* INSW(TSUMCROP - TSUMSBR, 0.,1.)

* RECOV, estimated below as a parameter, is the fraction of WCCR indicating favorable water content for the recovery from drought

* PUSHDORMRECTSUM is TSUM accumulation for the period in which PUSHDORMREC is ON

      PUSHDORMRECTSUM = INTGRL(ZERO, RPUSHDORMRECTSUM)

      RPUSHDORMRECTSUM = DTEFF*PUSHDORMREC - (PUSHDORMRECTSUM/DELT) * (1. - PUSHDORMREC) * (1. - PUSHREDIST)

* The dormancy ends when PUSHREDIST is ON (PUSHREDIST=1), triggering the start of DM distribution from storage roots to the
* leaves after a delay period of TSUM accumulation equivalent to DELREDIST (TSUM delay before DM distribution starts)
* PUSHREDIST is activated only when the condition for PUSHDORMREC is continuously favourable for the period of the delay.

      PUSHREDIST = INSW(PUSHDORMRECTSUM - DELREDIST,0.,1.)*(1. - PUSHREDISTEND)

* Duration of the DM distribution from storage roots to leaves (from PUSHREDIST=1,PUSHREDISTEND=0 to PUSHREDIST=0,PUSHREDISTEND=1)

      PUSHREDISTSUM = INTGRL(ZERO, RPUSHREDISTSUM)
      RPUSHREDISTSUM = DTEFF*PUSHREDIST - (PUSHREDISTSUM/DELT)*PUSHREDISTEND

* DM distribution from storage roots to leaves ends only when at least one of the following conditions are met:
* WSOREDISTFRAC (the fraction of DM redistributed) has reached a max DM redistribution fraction WSOREDISTFRACMAX, or
* WLVGREDISTCUMUL (the amount of new leaves produced) has reached WLVGNEWN (a min amount leaf DM that can be produced from
* the redistribution or  PUSHREDISTSUM (TSUM accumulation during the redistribution) has reached TSUMREDISTMAX (a max TSUM value
* for DM redistribution from storage roots). All these conditions can be met only after dormancy.

      PUSHREDISTEND = MAX(INSW(WSOREDISTFRAC-WSOREDISTFRACMAX,0.,1.), INSW(WLVGREDISTCUMUL-WLVGNEWN,0.,1.), ...
                    INSW(PUSHREDISTSUM-TSUMREDISTMAX,0.,1.)) * INSW(-PUSHREDISTSUM, 1.,0.)

      PUSHREDISTENDTSUM = INTGRL(ZERO, RPUSHREDISTENDTSUM)
      RPUSHREDISTENDTSUM = DTEFF*PUSHREDIST - (PUSHREDISTENDTSUM/DELT) * (1. - PUSHREDISTEND)

      WSOREDISTFRAC = (WSOREDISTCUMUL / NOTNUL(WSO))

      WSOREDISTCUMUL = INTGRL(ZERO, RREDISTSO)

* WSOREDISTCUMUL is the cumulative amount of storage roots DM produced after the distribution

* After a long period of drought and dormancy, part of storage roots DM is converted into energy
* for the production of new leaves. RREDISTSO is the rate of convertion of storage organs DM into new leaves


**** 6c. Dry matter redistribution after dormancy
** RREDISTLVG is the growth rate of new leaves with DM provided by storage roots: [g leaves DM/m2/day]
** SO2LV is the converting factor of the DM from storage roots to leaves : [g leaves DM/g storage roots DM]
** RREDISTSO is the rate of redistribution of DM from storage roots to leaves: [g storage roots DM/m2/day]
** WLVGREDISTCUMUL is the cumulative amount of new leaves DM produced after the distribution

     RREDISTLVG = SO2LV * RREDISTSO * (1. - DORMANCY)

     RREDISTSO = RRREDISTSO * WSO * PUSHREDIST - WSOREDISTCUMUL/DELT * INSW(-DORMTSUM, 1.,0.)

     WLVGREDISTCUMUL = INTGRL(ZERO, RREDISTLVG)

* The model assumes some energy or DM loss in producing new leaves from storage roots DM. This will be estimated
* using RREDISTMAINTLOSS

     RREDISTMAINTLOSS = (1. - SO2LV) * RREDISTSO

     REDISTMAINTLOSSCUMUL = INTGRL(ZERO, RREDISTMAINTLOSS)


***   8. Functions and parameters for cassava

*     Section 1
PARAM WCI   = 0.41
PARAM SLAI = 0.017; SLA_MAX = 0.03; NCUTTINGS = 1.5625; WCUTTINGUNIT = 14.
* Make sure to change SLAI as the first value of FRACSLATB x SLA_MAX whenever SLA_MAX is changed in sensitivity analysis.
* FRACSLATB is the fraction of SLA_MAX at different physiological times

*     Section 2
PARAM TBASE = 15.

*     Section 3
PARAM DOYPL = 113.

** New parameters
PARAM OPTEMERGTSUM = 180.; LAIEXPOEND = 0.75; TSUMLA_MIN = 180.; RGRL = 0.003; TSUMSOBULKINIT = 540.; RDRSHM = 0.09; RDRB = 0.09;
PARAM TSUMSBR = 780.; RDRWCUTTING = 0.017; LAICR = 3.5; LAI_MIN = 0.09; TSUMLLIFE = 1200.; FRACTLLFENHSH = 0.85; DELREDIST = 12.;
PARAM WSOREDISTFRACMAX = 0.05; WLVGNEWN = 10.; SO2LV = 0.8; RRREDISTSO = 0.01; FASTRANSLSO = 0.45; TSUMRROOTEND = 720.;
PARAM TSUMREDISTMAX = 144.; WCUTTINGMINPRO = 0.15; RECOV = 0.7; FST_CUTT = 0.03; FRT_CUTT = 0.03; FLV_CUTT = 0.07; FSO_CUTT = 0.

* When TSUMSOBULKINIT is changed, its value in the partitioning table should also be adjusted
* When TBASE is changed, RDRT and TTB should also be adjusted

FUNCTION FRACSLATB = 0.,0.57, 1440.,0.57, 2880.,0.65, 3864.,1., 4320.,1.

FUNCTION RDRT = -10.,0.02, 10.,0.02, 15.,0.03, 30.,0.06, 50.,0.06

FUNCTION TTB = -10.,0.0, 15.,0.0, 25.,1., 29.,1., 40.,0., 50.,0.0

*     Section 4
PARAM LUE_OPT = 1.5; K_ext = 0.67

*     Section 5
PARAM ROOTDM = 0.7; RRDMAX = 0.012

* ROOTDM = 2.6m (Connor et al., 1981, in: Van Heemst, 1988)

*     Partitioning tables for leaves (LV), stems (ST), storage organs (SO) and roots (RT):
**    (based on TSUMCROP, meaning TSUM after emergence)
FUNCTION FRTTB =     0.,0.11,   540.,0.100,  720.,0.094,  900.,0.01,  1488.,0.01, ...
      1980.,0.01, 2676.,0.01,  3864.,0.010, 4320.,0.010
FUNCTION FLVTB =     0.,0.71,   540.,0.515,  720.,0.393,  900.,0.24,  1488.,0.21, ...
      1980.,0.18, 2676.,0.13,  3864.,0.210, 4320.,0.210
FUNCTION FSTTB =     0.,0.18,   540.,0.385,  720.,0.393,  900.,0.26,  1488.,0.26, ...
      1980.,0.19, 2676.,0.29,  3864.,0.290, 4320.,0.290
FUNCTION FSOTB =     0.,0.00,   540.,0.000,  720.,0.120,  900.,0.49,  1488.,0.52, ...
      1980.,0.62, 2676.,0.57,  3864.,0.490, 4320.,0.490

**     Section 8
PARAM WCAD = 0.01; WCWP = 0.12; WCFC = 0.41; WCWET = 0.46; WCST = 0.52; TWCSD = 1.05; FRACRNINTC = 0.25
PARAM TRANCO = 8.; DRATE = 50.; IRRIGF = 0.

***********************************************************************
*Run 0: Sevekpota 2013
END

*Run 0: Sevekpota 2013
PARAM IRRIGF = 1.
END

STOP

* ---------------------------------------------------------------------*
*  SUBROUTINE GLA                                                      *
*  Purpose: This subroutine computes daily increase of leaf area index *
*           (ha leaf/ ha ground/ d)                                    *
* ---------------------------------------------------------------------*

      SUBROUTINE GLA(DTEFF,TSUMCROP,LAII,RGRL,DELT,SLA,LAI,
     $               GLV,TSUMLA_MIN,TRANRF,WC,WCWP,
     $               RWCUTTING,FLV,LAIEXPOEND,DORMANCY,
     $               GLAI)
      IMPLICIT REAL (A-Z)

*---- Growth during maturation stage:

** Accounting for dry matter redistribution from storage roots after dormancy

      GLAI = SLA * GLV * (1. - DORMANCY)

*---- Growth during juvenile stage:

** The model assumes the juvenile stage accounts for the period in which shoot and root growths depend on the reserves
** of the stem cutting and on photosynthesis for a duration of TSUMLA_MIN
      IF ((TSUMCROP.LT.TSUMLA_MIN).AND.(LAI.LT.LAIEXPOEND))
     $   GLAI = ((LAI * (EXP(RGRL * DTEFF * DELT) - 1.) / DELT) + ABS(RWCUTTING) * FLV * SLA) * TRANRF

*---- Growth at day of seedling emergence:
      IF ((TSUMCROP.GT.0.).AND.(LAI.EQ.0.).AND.(WC.GT.WCWP))
     $   GLAI = LAII / DELT

*---- Growth before seedling emergence:
      IF (TSUMCROP.EQ.0.)
     $   GLAI = 0.

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE PENMAN                                                   *
*  Purpose: Computation of the PENMAN EQUATION                         *
* ---------------------------------------------------------------------*

      SUBROUTINE PENMAN(DAVTMP,VP,DTR,LAI,WN,RNINTC,
     $                  PEVAP,PTRAN)
      IMPLICIT REAL (A-Z)

      DTRJM2 = DTR * 1.E6
      BOLTZM = 5.668E-8
      LHVAP  = 2.4E6
      PSYCH  = 0.067

      BBRAD  = BOLTZM * (DAVTMP+273.)**4 * 86400.
      SVP    = 0.611 * EXP(17.4 * DAVTMP / (DAVTMP + 239.))
      SLOPE  = 4158.6 * SVP / (DAVTMP + 239.)**2
      RLWN   = BBRAD * MAX(0.,0.55*(1.-VP/SVP))
      NRADS  = DTRJM2 * (1.-0.15) - RLWN
      NRADC  = DTRJM2 * (1.-0.25) - RLWN
      PENMRS = NRADS * SLOPE/(SLOPE+PSYCH)
      PENMRC = NRADC * SLOPE/(SLOPE+PSYCH)

      WDF    = 2.63 * (1.0 + 0.54 * WN)
      PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SLOPE+PSYCH)

      PEVAP  =     EXP(-0.5*LAI)  * (PENMRS + PENMD) / LHVAP
      PTRAN  = (1.-EXP(-0.5*LAI)) * (PENMRC + PENMD) / LHVAP
      PTRAN  = MAX( 0., PTRAN-0.5*RNINTC )

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE EVAPTR                                                   *
*  Purpose: To compute actual rates of evaporation and transpiration   *
*  WCCR was modified to ensure in order to never become                *
*  smaller than WCSD                                                   *
* ---------------------------------------------------------------------*

      SUBROUTINE EVAPTR(PEVAP,PTRAN,ROOTD,WA,WCAD,WCWP,WCFC,WCWET,WCST,
     $                  TRANCO,DELT,WCSD,
     $                  WCCR,EVAP,TRAN)
      IMPLICIT REAL (A-Z)

      WC   = 0.001 * WA   / ROOTD
      WAAD = 1000. * WCAD * ROOTD
      WAFC = 1000. * WCFC * ROOTD

      EVAP  = PEVAP * LIMIT( 0., 1., (WC-WCAD)/(WCFC-WCAD) )
         WCCR = WCWP + MAX(WCSD-WCWP, PTRAN/(PTRAN+TRANCO)*(WCFC-WCWP))

         IF (WC.GT.WCCR) THEN
             FR = LIMIT( 0., 1., (WCST-WC)/(WCST-WCWET) )
         ELSE
             FR = LIMIT( 0., 1., (WC-WCWP)/(WCCR-WCWP)  )
         ENDIF
      TRAN = PTRAN * FR

         AVAILF = MIN( 1., ((WA-WAAD)/DELT)/NOTNUL(EVAP+TRAN) )
      EVAP = EVAP * AVAILF
      TRAN = TRAN * AVAILF

      RETURN
      END

* ---------------------------------------------------------------------*
*  SUBROUTINE DRUNIR                                                   *
*  Purpose: To compute rates of drainage, runoff and irrigation        *
* ---------------------------------------------------------------------*

      SUBROUTINE DRUNIR(RAIN,RNINTC,EVAP,TRAN,IRRIGF,
     $                  DRATE,DELT,WA,ROOTD,WCFC,WCST,
     $                  DRAIN,RUNOFF,IRRIG)
      IMPLICIT REAL (A-Z)

      WC   = 0.001 * WA   / ROOTD
      WAFC = 1000. * WCFC * ROOTD
      WAST = 1000. * WCST * ROOTD

      DRAIN  = LIMIT( 0., DRATE, (WA-WAFC)/DELT +
     $               (RAIN - RNINTC - EVAP - TRAN)                  )

      RUNOFF =          MAX( 0., (WA-WAST)/DELT +
     $               (RAIN - RNINTC - EVAP - TRAN - DRAIN)          )

      IRRIG  = IRRIGF * MAX( 0., (WAFC-WA)/DELT -
     $               (RAIN - RNINTC - EVAP - TRAN - DRAIN - RUNOFF) )

      RETURN
      END

* ---------------------------------------------------------------------*
ENDJOB
