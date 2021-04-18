$title Affine Adjustable robust counterpart UC
$OnEmpty
$OffDigit

********************************************************************************
********************************CLI Configuation********************************
********************************************************************************
*
* usercuts: use User cuts algorith for N-1 constriants [default=0]
*   Note taht setting to 0 means that al contingencies will be considered
* renewable: turn on renewable generatos [default=1]
* past: past time periods of information to be considered [dafault=1]
* compactBal: use compact balance restricction, sum(gammas) = +-1 [default=1]
* deterministic: transform the model as a deterministic one [default=0]
* RC: transform the model as a robust counter part (non affine nor adjustable) [default=0]

$if not set usercuts $set usercuts 0
$if not set renewable $set renewable 1
$if not set past $set past 1
$if not set compactBal $set compactBal 1
$if not set deterministic $set deterministic 0
$if not set RC $set RC 0

$ifE %deterministic% == 1 $set compactBal 0
$ifE %RC% == 1 $set compactBal 0

scalar starttime, time_elapsed;
starttime = jnow;

********************************************************************************
**************************************sets**************************************
********************************************************************************

SETS
t        time periods                    /t1*t24/
l        lines                           /l1*l7/
i        bus                             /1*6/

phiHT    conventional generators         /j1*j3/
phiW     wind generators                 /q1/
phiS     PV generators                   /p1*p3/
;

alias(l, ll, k);
alias(t, tt);
alias(i, ii);
alias(phiHT, phiHHTT);

********************************************************************************
*************************input data and others**********************************
********************************************************************************

scalar idSlack idBarraReferencia /1/;
scalar lineCap capacidad de línea /1.1/;

TABLE genData(phiHT, *) conventional generator information
$include genData.inc
;

TABLE sData(phiS, *) solar data
$include sData.inc
;

TABLE lineData(l, *)
$include lineData.inc
;

TABLE d(t, i) demanda [% nivel de división de la carga]
$include d.inc
;

* source http://atlas.ideam.gov.co/visorAtlasVientos.html
table wOut(t, phiW) salida de generación eolica [MW]
$include wOut.inc
;

* source http://hrudnick.sitios.ing.uc.cl/alumno14/Hidroaysen/styled/index.html
table sOut(t, phiS) salida de generacion solar
$include sOut.inc
;

table gen_map(phiHT, i) lo que dice el nombre
   1  2  3  4  5  6
j1 1  0  0  0  0  0
j2 0  1  0  0  0  0
j3 0  0  0  0  0  1
;

table w_map(phiW, i) lo que dice el nombre
   1  2  3  4  5  6
q1 0  0  1  0  0  0
;

table s_map(phiS, i) lo que dice el nombre
   1  2  3  4  5  6
p1 0  0  0  1  0  0
p2 0  0  0  0  1  0
p3 0  0  1  0  0  0
;

scalar penalRenovable /300/;

wOut(t, phiW) = wOut(t, phiW)${%renewable% eq 1};
sOut(t, phiS) = sOut(t, phiS)${%renewable% eq 1};

********************************************************************************
**************************shift factors and related*****************************
********************************************************************************

Table PTDF(l, i) 'Matriz de PTDFS del caso base'
$include 6BUS_PTDF.inc
;

Table LODFS(l, ll) 'Matriz completa de LODFS de todo el Sistema'
$include 6BUS_LODF.inc
;

Table OTDF(ll, l, i) 'Matriz OTDF'
$include 6BUS_kptdf.inc
;

*otdf can be calculated as shown
*OTDF(k, l, i) = PTDF(l, i) + LODFS(l, k)*PTDF(k, i);

********************************************************************************
***********************User cuts parameter definitions**************************
********************************************************************************

parameters
PFijctg(t, l, k)  Power flow under contigency k
ContBin(t, l, k)  Auxiliar variable to activate user cuts
overLoad(t, l, k) Power overload of the line l
sumOverLoad       Sum of the power overload
sumcont           Sum of the ContBin parameter
tol               tolerance                                      /1e-6/
;

ContBin(t, l, k) = 1${%usercuts% eq 0};

********************************************************************************
*******************parameters for uncertainty data representation***************
********************************************************************************
parameters
eWindMax(t, phiW)        max wind error for each time period and wind generator
eSolarMax(t, phiS)       max solar error for time period and solar generator
eLoadMax(t, i)           max load error for each time period and bus
;

* 10% of max forecast deviation for renewable sources and 5% for load

eWindMax(t, phiW) = wOut(t, phiW)*0.1${%deterministic% eq 0};
eSolarMax(t, phiS) = sOut(t, phiS)*0.1${%deterministic% eq 0};
eLoadMax(t, i) = d(t, i)*0.05${%deterministic% eq 0};

********************************************************************************
*************************variables for aarc formulation*************************
********************************************************************************

variables
gammaWind(phiHT, t, phiW, tt) sensibilidad eolica
gammaSolar(phiHT, t, phiS, tt) sensibilidad solar
gammaLoad(phiHT, t, tt) sensilidad carga
;

** Non anticipativity conditions
gammaWind.fx(phiHT, t, phiW, tt)${ord(tt) gt ord(t)} = 0;
gammaSolar.fx(phiHT, t, phiS, tt)${ord(tt) gt ord(t)} = 0;
gammaLoad.fx(phiHT, t, tt)${ord(tt) gt ord(t)} = 0;

** limit of the past information that has been made certaint
gammaWind.fx(phiHT, t, phiW, tt)${abs(ord(t) - ord(tt)) gt %past%} = 0;
gammaSolar.fx(phiHT, t, phiS, tt)${abs(ord(t) - ord(tt)) gt %past%} = 0;
gammaLoad.fx(phiHT, t, tt)${abs(ord(t) - ord(tt)) gt %past%} = 0;


if(%deterministic% eq 1 or %RC% eq 1,
         gammaWind.fx(phiHT, t, phiW, tt) = 0;
         gammaSolar.fx(phiHT, t, phiS, tt) = 0;
         gammaLoad.fx(phiHT, t, tt) = 0;
);


********************************************************************************
******************************variable definitions******************************
********************************************************************************

Variables
zobj             Objective variable
pf(t, l)         power flow
pfk(t, l, ll)    power flow under contngencies
pnet(t, i)       net power per bus
;

binary variables
x(t, phiHT)      generator status variables
SU(t, phiHT)     startup variables
SD(t, phiHT)     shutdown variable
;

positive variables
p(t, phiHT)      active power of generatos
r(t, i)          load shedding (not used)
;

r.fx(t, i) = 0;

********************************************************************************
*************************variable definitions for aaro**************************
********************************************************************************

Variables

UbalWind(t, phiW, tt)            Auxiliar var for  power balance constraints
UbalSola(t, phiS, tt)            Auxiliar var for  power balance constraints
UbalLoad(t, tt)                  Auxiliar var for  power balance constraints

UabsLoad(l, t, i, tt)            Auxiliar var for power flow constraints
UabsWind(l, t, phiW, tt)         Auxiliar var for power flow constraints
UabsSola(l, t, phiS, tt)         Auxiliar var for power flow constraints

UabsLoadCont(l, t, i, tt, k)     Auxiliar var for n-1 power flow constraints
UabsWindCont(l, t, phiW, tt, k)  Auxiliar var for n-1 power flow constraints
UabsSolaCont(l, t, phiS, tt, k)  Auxiliar var for n-1 power flow constraints

UgenWind(phiHT, t, phiW, tt)     Auxiliar var for generation constraints
UgenSolar(phiHT, t, phiS, tt)    Auxiliar var for generation constraints
UgenLoad(phiHT, t, tt)           Auxiliar var for generation constraints

UfoAbsWind(phiW, tt)             Auxiliar var for obj. fun. constraints
UfoAbsSola(phiS, tt)             Auxiliar var for obj. fun. constraints
UfoAbsLoad(tt)                   Auxiliar var for obj. fun. constraints

URampWind(phiHT, t, phiW, tt)    Auxiliar var for ramping constraints
URampSolar(phiHT, t, phiS, tt)   Auxiliar var for ramping constraints
URampLoad(phiHT, t, tt)          Auxiliar var for ramping constraints
;

********************************************************************************
******************************Equation declaration******************************
********************************************************************************

Equations

costo                            objective function

maxGen(t, phiHT)                 Maximun generation
minGen(t, phiHT)                 Minimun generation

RampDo(t, phiHT)                 max ramp down
RampUp(t, phiHT)                 max ramp up

bal1(t, i)                       net power balance per bus
bal2(t)                          power balance

flujo(t, l)                      power flow
flujoMax(t, l)                   max power flow
flujoMin(t, l)                   min power flow

flujoCont(t, l, k)               power flow under contingency
flujoMaxCont(t, l, k)            max power flow under contingency
flujoMinCont(t, l, k)            min power flow under contingency

status1(t, phiHT)                on-off logic
status2(t, phiHT)                non simultaneous starup and shutdown
min_updown_1(t, phiHT)           min up-down time at t=0
min_updown_2(t,phiHT)            min up time
min_updown_3(t,phiHT)            min down time

********************************************************************************
**************************aditional equations for aarc**************************
********************************************************************************

UabsWindMaxCont(l, t, phiW, tt, k) Aux. equation for N-1 power flow limits
UabsWindMinCont(l, t, phiW, tt, k) Aux. equation for N-1 power flow limits
UabsSolaMaxCont(l, t, phiS, tt, k) Aux. equation for N-1 power flow limits
UabsSolaMinCont(l, t, phiS, tt, k) Aux. equation for N-1 power flow limits
UabsLoadMaxCont(l, t, i, tt, k)    Aux. equation for N-1 power flow limits
UabsLoadMinCont(l, t, i, tt, k)    Aux. equation for N-1 power flow limits

UbalWindMax(t, phiW, tt)           Aux. equation for power balance
UbalWindMin(t, phiW, tt)           Aux. equation for power balance
UbalSolaMax(t, phiS, tt)           Aux. equation for power balance
UbalSolaMin(t, phiS, tt)           Aux. equation for power balance
UbalLoadMax(t, tt)                 Aux. equation for power balance
UbalLoadMin(t, tt)                 Aux. equation for power balance

UabsWindMax(l, t, phiW, tt)        Aux. equation for power flow limits
UabsWindMin(l, t, phiW, tt)        Aux. equation for power flow limits
UabsSolaMax(l, t, phiS, tt)        Aux. equation for power flow limits
UabsSolaMin(l, t, phiS, tt)        Aux. equation for power flow limits
UabsLoadMax(l, t, i, tt)           Aux. equation for power flow limits
UabsLoadMin(l, t, i, tt)           Aux. equation for power flow limits

UgenAbsWindMax(phiHT, t, phiW, tt) Aux. equation for power generator limits
UgenAbsWindMin(phiHT, t, phiW, tt) Aux. equation for power generator limits
UgenAbsSolaMax(phiHT, t, phiS, tt) Aux. equation for power generator limits
UgenAbsSolaMin(phiHT, t, phiS, tt) Aux. equation for power generator limits
UgenAbsLoadMax(phiHT, t, tt)       Aux. equation for power generator limits
UgenAbsLoadMin(phiHT, t, tt)       Aux. equation for power generator limits

UfoWindMax(phiW, tt)               Aux. equation for objective function
UfoWindMin(phiW, tt)               Aux. equation for objective function
UfoSolaMax(phiS, tt)               Aux. equation for objective function
UfoSolaMin(phiS, tt)               Aux. equation for objective function
UfoLoadMax(tt)                     Aux. equation for objective function
UfoLoadMin(tt)                     Aux. equation for objective function

URampWindUp(phiHT, t, phiW, tt)    Aux. equation for ramping constraints
URampWindDn(phiHT, t, phiW, tt)    Aux. equation for ramping constraints
URampSolarUp(phiHT, t, phiS, tt)   Aux. equation for ramping constraints
URampSolarDn(phiHT, t, phiS, tt)   Aux. equation for ramping constraints
URampLoadUp(phiHT, t, tt)          Aux. equation for ramping constraints
URampLoadDn(phiHT, t, tt)          Aux. equation for ramping constraints

** alternative power balance
newbal(t)                         compact form of the uncertainty power balance
newGW(t, phiW, tt)                compact form of the uncertainty power balance
newGS(t, phiS, tt)                compact form of the uncertainty power balance
newGL(t, tt)                      compact form of the uncertainty power balance

;

********************************************************************************
*******************************Equation definition******************************
********************************************************************************

***************************************MIP**************************************

status1(t, phiHT)..
         x(t, phiHT) - x(t-1, phiHT) =e= SU(t, phiHT) - SD(t, phiHT)
;

status2(t, phiHT)..
         SU(t, phiHT) + SD(t, phiHT) =l= 1
;

min_updown_1(t, phiHT)$(genData(phiHT, 'MinUp') + genData(phiHT, 'MinDn') gt 0 and ord(t) le genData(phiHT, 'MinUp') + genData(phiHT, 'MinDn'))..
         x(t, phiHT) =e= 1${genData(phiHT, 't0')};


min_updown_2(t,phiHT)..
         sum(tt$(ord(tt) ge ord(t)-genData(phiHT, 'MinUp')+1 and ord(tt) le ord(t)), SU(t, phiHT)) =l= x(t,phiHT);

min_updown_3(t,phiHT)..
         sum(tt$(ord(tt) ge ord(t)-genData(phiHT, 'MinDn')+1 and ord(tt) le ord(t)), SD(t, phiHT)) =l= 1-x(t,phiHT);

*******************************objective function*******************************

costo..
         zobj
         =g=
         sum((phiHT, t),
                 p(t, phiHT)*genData(phiHT, 'Cost1')*genData(phiHT, 'fp')
                 + x(t, phiHT)*genData(phiHT, 'Cost0')*genData(phiHT, 'fp')
         )
         + sum((t, phiHt), SU(t, phiHT)*genData(phiHT, 'Csu') + SD(t, phiHT)*genData(phiHT, 'Csd'))
** AARC additional term for objective function
         + sum(tt,
                 + sum(phiW, UfoAbsWind(phiW, tt)*eWindMax(tt, phiW))
                 + sum(phiS, UfoAbsSola(phiS, tt)*eSolarMax(tt, phiS))
                 + sum(i,    UfoAbsLoad(tt)*eLoadMax(tt, i))
         )
;

** Absolute values needed on the computation of the objective function
UfoWindMax(phiW, tt)..
         UfoAbsWind(phiW, tt)
         =g=
         sum((phiHT, t)${ord(tt) le ord(t)},
                 gammaWind(phiHT, t, phiW, tt)*genData(phiHT, 'Cost1')*genData(phiHT, 'fp')
         )
;

UfoWindMin(phiW, tt)..
         -UfoAbsWind(phiW, tt)
         =l=
         sum((phiHT, t)${ord(tt) le ord(t)},
                  gammaWind(phiHT, t, phiW, tt)*genData(phiHT, 'Cost1')*genData(phiHT, 'fp')
         )
;

UfoSolaMax(phiS, tt)..
         UfoAbsSola(phiS, tt)
         =g=
         sum((phiHT, t)${ord(tt) le ord(t)},
                  gammaSolar(phiHT, t, phiS, tt)*genData(phiHT, 'Cost1')*genData(phiHT, 'fp')
         )
;

UfoSolaMin(phiS, tt)..
         -UfoAbsSola(phiS, tt)
         =l=
         sum((phiHT, t)${ord(tt) le ord(t)},
                 gammaSolar(phiHT, t, phiS, tt)*genData(phiHT, 'Costo1')*genData(phiHT, 'fp')
         )
;

UfoLoadMax(tt)..
         UfoAbsLoad(tt)
         =g=
         sum((phiHT, t)${ord(tt) le ord(t)},
                  gammaLoad(phiHT, t, tt)*genData(phiHT, 'Costo1')*genData(phiHT, 'fp')
         )
;

UfoLoadMin(tt)..
         -UfoAbsLoad(tt)
         =l=
         sum((phiHT, t)${ord(tt) le ord(t)},
                  gammaLoad(phiHT, t, tt)*genData(phiHT, 'Costo1')*genData(phiHT, 'fp')
         )
;

**********************************power balance*********************************

bal1(t, i)..
         pnet(t, i)
         =e=
                 + sum[phiHT, gen_map(phiHT, i) * p(t, phiHT)]
                 + sum[phiW , w_map(phiW, i) * wOut(t, phiW)]
                 + sum[phiS , s_map(phiS, i) * sOut(t, phiS)]
                 - d(t, i)
                 + R(t, i)
;

bal2(t)..
         sum(i, pnet(t, i))
         - sum(tt${ord(tt) le ord(t)},
                 + sum(phiW, UbalWind(t, phiW, tt)*eWindMax(tt, phiW))
                 + sum(phiS, UbalSola(t, phiS, tt)*eSolarMax(tt, phiS))
                 + sum(i, UbalLoad(t, tt)*eLoadMax(tt, i))
         )
         =g= 0
;

UbalWindMax(t, phiW, tt)${ord(tt) le ord(t)}..
         sum(phiHT, gammaWind(phiHT, t, phiW, tt)) + 1${ord(tt) eq ord(t)} =l= UbalWind(t, phiW, tt)
;

UbalWindMin(t, phiW, tt)${ord(tt) le ord(t)}..
         -UbalWind(t, phiW, tt) =l= sum(phiHT, gammaWind(phiHT, t, phiW, tt)) + 1${ord(tt) eq ord(t)}
;

UbalSolaMax(t, phiS, tt)${ord(tt) le ord(t)}..
         sum(phiHT, gammaSolar(phiHT, t, phiS, tt)) + 1${ord(tt) eq ord(t)} =l= UbalSola(t, phiS, tt)
;

UbalSolaMin(t, phiS, tt)${ord(tt) le ord(t)}..
         -UbalSola(t, phiS, tt) =l= sum(phiHT, gammaSolar(phiHT, t, phiS, tt)) + 1${ord(tt) eq ord(t)}
;

UbalLoadMax(t, tt)${ord(tt) le ord(t)}..
         sum(phiHT, gammaLoad(phiHT, t, tt)) - 1${ord(tt) eq ord(t)} =l= UbalLoad(t, tt)
;

UbalLoadMin(t, tt)${ord(tt) le ord(t)}..
         -UbalLoad(t, tt) =l= sum(phiHT, gammaLoad(phiHT, t, tt)) - 1${ord(tt) eq ord(t)}
;

***************compact for of the power balance under uncertainty***************

newbal(t)..
         + sum(phiHHTT, p(t, phiHHTT))
         + sum(phiW, wOut(t, phiW))
         + sum(phiS, sOut(t, phiS)*sData(phiS, 'Pmax'))
         - sum(i, d(t, i) - R(t, i))
         =e=
         0
;

newGW(t, phiW, tt)${ord(tt) le ord(t)}..
         sum(phiHT, gammaWind(phiHT, t, phiW, tt)) =e= - 1${ord(tt) eq ord(t)}
;

newGS(t, phiS, tt)${ord(tt) le ord(t)}..
         sum(phiHT, gammaSolar(phiHT, t, phiS, tt)) =e= - 1${ord(tt) eq ord(t)}
;

newGL(t, tt)${ord(tt) le ord(t)}..
         sum(phiHT, gammaLoad(phiHT, t, tt)) =e= + 1${ord(tt) eq ord(t)}
;

***********************************power flow***********************************

flujo(t, l)..
         pf(t, l)
         =e=
         sum(i, PTDF(l, i)*pnet(t, i))
;

flujoMax(t, l)..
         pf(t, l)
         =l=
         lineData(l, 'limit')*lineCap
         + sum(tt${ord(tt) le ord(t)},
                  - sum(phiS, UabsSola(l, t, phiS, tt)*eSolarMax(tt, phiS))
                  - sum(phiW, UabsWind(l, t, phiW, tt)*eWindMax(tt, phiW))
                  - sum(i, UabsLoad(l, t, i, tt)*eLoadMax(tt, i))
         )
;

flujoMin(t, l)..
         - pf(t, l)
         =l=
         lineData(l, 'limit')*lineCap
         + sum(tt${ord(tt) le ord(t)},
                  - sum(phiS, UabsSola(l, t, phiS, tt)*eSolarMax(tt, phiS))
                  - sum(phiW, UabsWind(l, t, phiW, tt)*eWindMax(tt, phiW))
                  - sum(i, UabsLoad(l, t, i, tt)*eLoadMax(tt, i))
         )
;

UabsWindMax(l, t, phiW, tt)${ord(tt) le ord(t)}..
         - UabsWind(l, t, phiW, tt)
         =l=
         sum(i,
                 PTDF(l, i)
                 * (
                         + w_map(phiW, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                 gen_map(phiHT, i)*gammaWind(phiHT, t, phiW, tt)
                         )
                 )
         )
;

UabsWindMin(l, t, phiW, tt)${ord(tt) le ord(t)}..
         sum(i,
                 PTDF(l, i)
                 * (
                         + w_map(phiW, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                  gen_map(phiHT, i)*gammaWind(phiHT, t, phiW, tt)
                         )
                 )
         )
         =l=
         UabsWind(l, t, phiW, tt)
;

UabsSolaMax(l, t, phiS, tt)${ord(tt) le ord(t)}..
         - UabsSola(l, t, phis, tt)
         =l=
         sum(i,
                 PTDF(l, i)
                 * (
                         + s_map(phiS, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                 gen_map(phiHT, i)*gammaSolar(phiHT, t, phiS, tt)
                         )
                 )
         )
;

UabsSolaMin(l, t, phiS, tt)${ord(tt) le ord(t)}..
         sum(i,
                 PTDF(l, i)
                 * (
                         + s_map(phiS, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                 gen_map(phiHT, i)*gammaSolar(phiHT, t, phiS, tt)
                         )
                 )
         )
         =l=
         UabsSola(l, t, phiS, tt)
;

UabsLoadMax(l, t, i, tt)${ord(tt) le ord(t)}..
         - UabsLoad(l, t, i, tt)
         =l=
         sum(ii,
                 PTDF(l, ii)
                 *(
                         sum(phiHT,
                                 gen_map(phiHT, ii)*gammaLoad(phiHT, t, tt)
                         )
                         - 1${ord(t) eq ord(tt) and ord(i) eq ord(ii)}
                  )
         )
;

UabsLoadMin(l, t, i, tt)${ord(tt) le ord(t)}..
         sum(ii,
                 PTDF(l, ii)
                 *(
                         sum(phiHT,
                                 gen_map(phiHT, ii)*gammaLoad(phiHT, t, tt)
                         )
                         - 1${ord(t) eq ord(tt) and ord(i) eq ord(ii)}
                  )
         )
         =l=
         UabsLoad(l, t, i, tt)
;

**************************power flow under contingency**************************

flujoCont(t, l, k)${ContBin(t, l, k) ne 0}..
         pfk(t, l, k)
         =e=
         sum(i, OTDF(k, l, i)*pnet(t, i))
;

flujoMaxCont(t, l, k)${ContBin(t, l, k) ne 0}..
         pfk(t, l, k)
         =l=
         lineData(l, 'limit')*lineCap
         + sum(tt${ord(tt) le ord(t)},
                  - sum(phiS, UabsSolaCont(l, t, phiS, tt, k)*eSolarMax(tt, phiS))
                  - sum(phiW, UabsWindCont(l, t, phiW, tt, k)*eWindMax(tt, phiW))
                  - sum(i, UabsLoadCont(l, t, i, tt, k)*eLoadMax(tt, i))
         )
;

flujoMinCont(t, l, k)${ContBin(t, l, k) ne 0}..
         - pfk(t, l, k)
         =l=
         lineData(l, 'limit')*lineCap
         + sum(tt${ord(tt) le ord(t)},
                  - sum(phiS, UabsSolaCont(l, t, phiS, tt, k)*eSolarMax(tt, phiS))
                  - sum(phiW, UabsWindCont(l, t, phiW, tt, k)*eWindMax(tt, phiW))
                  - sum(i, UabsLoadCont(l, t, i, tt, k)*eLoadMax(tt, i))
         )
;

UabsWindMaxCont(l, t, phiW, tt, k)${ord(tt) le ord(t) and ContBin(t, l, k) ne 0}..
         - UabsWindCont(l, t, phiW, tt, k)
         =l=
         sum(i,
                 OTDF(k, l, i)
                 * (
                         + w_map(phiW, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                 gen_map(phiHT, i)*gammaWind(phiHT, t, phiW, tt)
                         )
                 )
         )
;

UabsWindMinCont(l, t, phiW, tt, k)${ord(tt) le ord(t) and ContBin(t, l, k) ne 0}..
         sum(i,
                 OTDF(k, l, i)
                 * (
                         + w_map(phiW, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                  gen_map(phiHT, i)*gammaWind(phiHT, t, phiW, tt)
                         )
                 )
         )
         =l=
         UabsWindCont(l, t, phiW, tt, k)
;

UabsSolaMaxCont(l, t, phiS, tt, k)${ord(tt) le ord(t) and ContBin(t, l, k) ne 0}..
         - UabsSolaCont(l, t, phis, tt, k)
         =l=
         sum(i,
                 OTDF(k, l, i)
                 * (
                         + s_map(phiS, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                 gen_map(phiHT, i)*gammaSolar(phiHT, t, phiS, tt)
                         )
                 )
         )
;

UabsSolaMinCont(l, t, phiS, tt, k)${ord(tt) le ord(t) and ContBin(t, l, k) ne 0}..
         sum(i,
                 OTDF(k, l, i)
                 * (
                         + s_map(phiS, i)${ord(t) eq ord(tt)}
                         + sum(phiHT,
                                 gen_map(phiHT, i)*gammaSolar(phiHT, t, phiS, tt)
                         )
                 )
         )
         =l=
         UabsSolaCont(l, t, phiS, tt, k)
;

UabsLoadMaxCont(l, t, i, tt, k)${ord(tt) le ord(t) and ContBin(t, l, k) ne 0}..
         - UabsLoadCont(l, t, i, tt, k)
         =l=
         sum(ii,
                 OTDF(k, l, ii)
                 *(
                         sum(phiHT,
                                 gen_map(phiHT, ii)*gammaLoad(phiHT, t, tt)
                         )
                         - 1${ord(t) eq ord(tt) and ord(i) eq ord(ii)}
                  )
         )
;

UabsLoadMinCont(l, t, i, tt, k)${ord(tt) le ord(t) and ContBin(t, l, k) ne 0}..
         sum(ii,
                 OTDF(k, l, ii)
                 *(
                         sum(phiHT,
                                 gen_map(phiHT, ii)*gammaLoad(phiHT, t, tt)
                         )
                         - 1${ord(t) eq ord(tt) and ord(i) eq ord(ii)}
                  )
         )
         =l=
         UabsLoadCont(l, t, i, tt, k)
;

*****************************power generation limits****************************

maxGen(t, phiHT)..
         p(t, phiHT)
         + sum(tt${ord(tt) le ord(t)},
                 + sum(phiW, UgenWind(phiHT, t, phiW, tt)*eWindMax(tt, phiW))
                 + sum(phiS, UgenSolar(phiHT, t, phiS, tt)*eSolarMax(tt, phiS))
                 + sum(i,    UgenLoad(phiHT, t, tt)*eLoadMax(tt, i))
         )
         =l=
         genData(phiHT, 'Pmax')*x(t, phiHT)
;

minGen(t, phiHT)..
         -p(t, phiHT)
         + sum(tt${ord(tt) le ord(t)},
                 + sum(phiW, UgenWind(phiHT, t, phiW, tt)*eWindMax(tt, phiW))
                 + sum(phiS, UgenSolar(phiHT, t, phiS, tt)*eSolarMax(tt, phiS))
                 + sum(i,    UgenLoad(phiHT, t, tt)*eLoadMax(tt, i))
         )
         =l=
         -genData(phiHT, 'Pmin')*x(t, phiHT)
;

UgenAbsWindMax(phiHT, t, phiW, tt)${ord(tt) le ord(t)}..
         -UgenWind(phiHT, t, phiW, tt) =l= gammaWind(phiHT, t, phiW, tt)
;

UgenAbsWindMin(phiHT, t, phiW, tt)${ord(tt) le ord(t)}..
         gammaWind(phiHT, t, phiW, tt) =l= UgenWind(phiHT, t, phiW, tt)
;

UgenAbsSolaMax(phiHT, t, phiS, tt)${ord(tt) le ord(t)}..
         -UgenSolar(phiHT, t, phiS, tt) =l= gammaSolar(phiHT, t, phiS, tt)
;

UgenAbsSolaMin(phiHT, t, phiS, tt)${ord(tt) le ord(t)}..
         gammaSolar(phiHT, t, phiS, tt) =l= UgenSolar(phiHT, t, phiS, tt)
;

UgenAbsLoadMax(phiHT, t, tt)${ord(tt) le ord(t)}..
         -UgenLoad(phiHT, t, tt) =l= gammaLoad(phiHT, t, tt)
;

UgenAbsLoadMin(phiHT, t, tt)${ord(tt) le ord(t)}..
         gammaLoad(phiHT, t, tt) =l= UgenLoad(phiHT, t, tt)
;

*******************************ramping constraints******************************

RampDo(t, phiHT)${ord(t) gt 1}..
         - genData(phiHT, 'Rdo')
         =l=
         p(t, phiHT) - p(t-1, phiHT)
         - sum(tt${ord(tt) le ord(t)},
                 + sum(phiW, URampWind(phiHT, t, phiW, tt)*eWindMax(tt, phiW))
                 + sum(phiS, URampSolar(phiHT, t, phiS, tt)*eSolarMax(tt, phiS))
                 + sum(i,    URampLoad(phiHT, t, tt)*eLoadMax(tt, i))
         )
;

RampUp(t, phiHT)${ord(t) gt 1}..
         p(t, phiHT) - p(t-1, phiHT)
         + sum(tt${ord(tt) le ord(t)},
                 + sum(phiW, URampWind(phiHT, t, phiW, tt)*eWindMax(tt, phiW))
                 + sum(phiS, URampSolar(phiHT, t, phiS, tt)*eSolarMax(tt, phiS))
                 + sum(i,    URampLoad(phiHT, t, tt)*eLoadMax(tt, i))
         )
         =l=
         genData(phiHT, 'Rup')
;

URampWindUp(phiHT, t, phiW, tt)${ord(tt) le ord(t)}..
         -URampWind(phiHT, t, phiW, tt)
         =l=
         + gammaWind(phiHT, t, phiW, tt)
         - gammaWind(phiHT, t-1, phiW, tt)${ord(tt) lt ord(t)}
;

URampWindDn(phiHT, t, phiW, tt)${ord(tt) le ord(t)}..
         + gammaWind(phiHT, t, phiW, tt)
         - gammaWind(phiHT, t-1, phiW, tt)${ord(tt) lt ord(t)}
         =l=
         URampWind(phiHT, t, phiW, tt)
;

URampSolarUp(phiHT, t, phiS, tt)${ord(tt) le ord(t)}..
         -URampSolar(phiHT, t, phiS, tt)
         =l=
         + gammaSolar(phiHT, t, phiS, tt)
         - gammaSolar(phiHT, t-1, phiS, tt)${ord(tt) lt ord(t)}
;

URampSolarDn(phiHT, t, phiS, tt)${ord(tt) le ord(t)}..
         + gammaSolar(phiHT, t, phiS, tt)
         - gammaSolar(phiHT, t-1, phiS, tt)${ord(tt) lt ord(t)}
         =l=
         URampSolar(phiHT, t, phiS, tt)
;

URampLoadUp(phiHT, t, tt)${ord(tt) le ord(t)}..
         -URampLoad(phiHT, t, tt)
         =l=
         + gammaLoad(phiHT, t, tt)
         - gammaLoad(phiHT, t-1, tt)${ord(tt) lt ord(t)}
;

URampLoadDn(phiHT, t, tt)${ord(tt) le ord(t)}..
         + gammaLoad(phiHT, t, tt)
         - gammaLoad(phiHT, t-1, tt)${ord(tt) lt ord(t)}
         =l=
         URampLoad(phiHT, t, tt)
;

********************************************************************************
* Model definitions
********************************************************************************

model UC /status1 ,status2,
*         min_updown_1,
         min_updown_2, min_updown_3/;
model costoTot /costo, UfoWindMax, UfoWindMin, UfoSolaMax, UfoSolaMin, UfoLoadMax, UfoLoadMin/;
model flujosNorm /flujo, flujoMax, flujoMin, UabsWindMax, UabsWindMin, UabsSolaMax, UabsSolaMin, UabsLoadMax, UabsLoadMin/;
model Continge /flujoCont, flujoMaxCont, flujoMinCont, UabsWindMaxCont, UabsWindMinCont, UabsSolaMaxCont, UabsSolaMinCont, UabsLoadMaxCont, UabsLoadMinCont/;
model rampConst /RampDo, RampUp, URampWindUp, URampWindDn, URampSolarUp, URampSolarDn, URampLoadUp, URampLoadDn/;
model limGen /maxGen, minGen, UgenAbsWindMax, UgenAbsWindMin, UgenAbsLoadMax, UgenAbsSolaMin, UgenAbsSolaMax, UgenAbsLoadMin/;

model balAfin /bal1, bal2, UbalWindMax, UbalWindMin, UbalSolaMax, UbalSolaMin, UbalLoadMax, UbalLoadMin/;
model balGammas /bal1, newbal, newGW, newGS, newGL/;

model basemodel /

* Base
         costoTot
         limGen
*         balGammas
*         balAfin
* Other restrictions
         UC
         rampConst
         flujosNorm
         Continge
         /;

model robust /basemodel, balGammas/;
model nocompactrobust /basemodel, balAfin/;

*********************************Configuration**********************************

* Configuration
option optcr = 1e-7;
option solveopt=replace;

option reslim = 86400;
option Savepoint=1;
robust.optfile = 1;
nocompactrobust.optfile=1;

*ep.solprint=1
* Turn off the listing of the input file *
$offlisting
* Turn off the listing and cross-reference of the symbols used
*$offsymxref offsymlist

option solveopt=replace;
option solvelink=2;

option
    limrow = 0
    limcol = 0
*    solprint = off
*    sysout = off
;
file opt cplex option file /cplex.opt/;
put opt;
put 'threads 0'/;
put 'miptrace mipAARC.csv'/;
putclose;

********************************************************************************
*******************************user cuts algorithm******************************
********************************************************************************

parameters
adjust(t, l, k)                  pesimistic total power flow addition
absWind(l, t, phiW, tt, k)       pesimistic power flow addition (wind)
absSola(l, t, phiS, tt, k)       pesimistic power flow addition (solar)
absLoad(l, t, i, tt, k)          pesimistic power flow addition (load)
;

repeat(

         overLoad(t, l, k) = 0;
         if(%compactBal% eq 1,
                 solve robust us mip min zobj;
         else
                 solve nocompactrobust us mip min zobj;
         );

**** pesimistic power flow calculation
         absWind(l, t, phiW, tt, k)${ord(tt) le ord(t)} =
         abs(
                 sum(i,
                         OTDF(k, l, i)
                         * (
                                 + w_map(phiW, i)${ord(t) eq ord(tt)}
                                 + sum(phiHT,
                                         gen_map(phiHT, i)*gammaWind.l(phiHT, t, phiW, tt)
                                 )
                         )
                 )
         )
         ;

         absSola(l, t, phiS, tt, k)${ord(tt) le ord(t)} =
         abs(
                 sum(i,
                         OTDF(k, l, i)
                         * (
                                 + s_map(phiS, i)${ord(t) eq ord(tt)}
                                 + sum(phiHT,
                                         gen_map(phiHT, i)*gammaSolar.l(phiHT, t, phiS, tt)
                                 )
                         )
                 )
         )
         ;

         absLoad(l, t, i, tt, k)${ord(tt) le ord(t)} =
         abs(
                 sum(ii,
                         OTDF(k, l, ii)
                         *(
                                 sum(phiHT,
                                         gen_map(phiHT, ii)*gammaLoad.l(phiHT, t, tt)
                                 )
                                 - 1${ord(t) eq ord(tt) and ord(i) eq ord(ii)}
                         )
                 )
         )
         ;

         adjust(t, l, k) =
         sum(tt${ord(tt) le ord(t)},
                  + sum(phiS, absSola(l, t, phiS, tt, k)*eSolarMax(tt, phiS))
                  + sum(phiW, absWind(l, t, phiW, tt, k)*eWindMax(tt, phiW))
                  + sum(i,    absLoad(l, t, i, tt, k)*eLoadMax(tt, i))
         )
         ;

         PFijctg(t, l, k) = abs[sum(i, OTDF(k, l, i)*pnet.l(t, i))] + adjust(t, l, k);

         ContBin(t, l, k)${PFijctg(t, l, k) gt lineData(l, 'limit')*lineCap} = 1;
         sumcont = sum((t, l, k), ContBin(t, l, k));

         overLoad(t, l, k)${PFijctg(t, l, k) gt lineData(l, 'limit')*lineCap} = PFijctg(t, l, k) - lineData(l, 'limit')*lineCap;
         sumOverLoad = sum((t, l, k) ,overLoad(t, l, k));

until {sumOverLoad} lt tol or %usercuts% eq 0);

********************************************************************************

parameters
solveTime
solverTime
foval
ggg(t, phiHT)
OutputGap;


foval = zobj.l;
solveTime = robust.etSolve;
solverTime = robust.etSolver;
ggg(t, phiHT) = eps + p.l(t, phiHT);
time_elapsed = (jnow - starttime)*24*3600;
OutputGap = (abs(robust.objest - robust.objval)) / max(abs(robust.objest), abs(robust.objval));

execute_unload 'SIM_past%past%.gdx';
execute_unload 'ans.gdx' foval ggg solveTime SolverTime time_elapsed OutputGap;
