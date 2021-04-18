$title deterministic formulation
$OnEmpty
$OffDigit

********************************************************************************
********************************CLI Configuation********************************
********************************************************************************
*
* usercuts: use User cuts algorith for N-1 constriants [default=0]
*   Note taht setting to 0 means that al contingencies will be considered
* renewable: turn on renewable generatos [default=1]

$if not set usercuts $set usercuts 0
$if not set renewable $set renewable 1

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
         sum(i, pnet(t, i)) =g= 0
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
;

flujoMin(t, l)..
         - pf(t, l)
         =l=
         lineData(l, 'limit')*lineCap
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
;

flujoMinCont(t, l, k)${ContBin(t, l, k) ne 0}..
         - pfk(t, l, k)
         =l=
         lineData(l, 'limit')*lineCap
;

*****************************power generation limits****************************

maxGen(t, phiHT)..
         p(t, phiHT)
         =l=
         genData(phiHT, 'Pmax')*x(t, phiHT)
;

minGen(t, phiHT)..
         -p(t, phiHT)
         =l=
         -genData(phiHT, 'Pmin')*x(t, phiHT)
;

*******************************ramping constraints******************************

RampDo(t, phiHT)${ord(t) gt 1}..
         - genData(phiHT, 'Rdo')
         =l=
         p(t, phiHT) - p(t-1, phiHT)
;

RampUp(t, phiHT)${ord(t) gt 1}..
         p(t, phiHT) - p(t-1, phiHT)
         =l=
         genData(phiHT, 'Rup')
;

********************************************************************************
* Model definitions
********************************************************************************

model UC /status1 ,status2,
*         min_updown_1,
         min_updown_2, min_updown_3/;
model costoTot /costo/;
model flujosNorm /flujo, flujoMax, flujoMin/;
model Continge /flujoCont, flujoMaxCont, flujoMinCont/;
model rampConst /RampDo, RampUp/;
model limGen /maxGen, minGen/;
model balance /bal1, bal2/;

model deterministic /

* Base
         costoTot
         limGen
         balance
         UC
         rampConst
         flujosNorm
         Continge
         /;

*********************************Configuration**********************************

* Configuration
option optcr = 1e-7;
option solveopt=replace;

option reslim = 86400;
option Savepoint=1;
deterministic.optfile = 1;

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
put 'miptrace mipDet.csv'/;
putclose;

********************************************************************************
*******************************user cuts algorithm******************************
********************************************************************************

repeat(

         overLoad(t, l, k) = 0;

         solve deterministic us mip min zobj;

         PFijctg(t, l, k) = abs[sum(i, OTDF(k, l, i)*pnet.l(t, i))];

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
solveTime = deterministic.etSolve;
solverTime = deterministic.etSolver;
ggg(t, phiHT) = eps + p.l(t, phiHT);
time_elapsed = (jnow - starttime)*24*3600;
OutputGap = (abs(deterministic.objest - deterministic.objval)) / max(abs(deterministic.objest), abs(deterministic.objval));

execute_unload 'deterministic.gdx';
execute_unload 'ans.gdx' foval ggg solveTime SolverTime time_elapsed OutputGap;
