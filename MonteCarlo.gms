$Title Montecarlo simulation
$OnEmpty
$OffDigit

********************************************************************************
********************************CLI Configuation********************************
********************************************************************************
*
* past: past time periods of information to be considered [dafault=1]

$if not set past $set past 1

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
******************************parameter definitios******************************
********************************************************************************

parameters
idSlack
lineCap
genData(phiHT, *)
gen_map(phiHT, i)
wOut(t, phiW)
sOut(t, phiS)
s_map(phiS, i)
w_map(phiW, i)
sData(phiS, *)
lineData(l, *)
d(t, i)
PTDF(l, i)

gslack(t)
costoSlack
ggg(t, phiHT)
Pneta(t, i)
SumNeta(t)
pga(t, phiHT)
costot
gammaSlackWind(t, phiW, tt)
gammaSlackSola(t, phiS, tt)
gammaSlackLoad(t, tt)

R(t, i)
p(t, phiHT)
x(t, phiHT)

gammaWind(phiHT, t, phiW, tt)
gammaSolar(phiHT, t, phiS, tt)
gammaLoad(phiHT, t, tt)

eWindMax(t, phiW)
eSolarMax(t, phiS)
eLoadMax(t, ii)

otdf(k, l, i)
;

** load the info from the gdx file

$gdxin 'SIM_past%past%.gdx'
$load idSlack lineCap genData gen_map wOut w_map sOut sData s_map lineData d PTDF R=R.l p=p.l x=x.l gammaWind = gammaWind.l gammaSolar = gammaSolar.l gammaLoad = gammaLoad.l eWindMax eSolarMax eLoadMax otdf
$gdxin

********************************************************************************
******************************MonteCarlo simulation*****************************
********************************************************************************

** random seed
execseed = 8347294;

sets
sim      MC simulations /mc1*mc100/
scen     scenarios set /i1*i300/;

scalars
bound    max devition on the simulations /0.0/
tol      tolerance for float point operations /1e-6/
;

parameters
pgA(t, phiHT)            Affine policy (a.k.a. linear or affine decision rule)
powBal(t)                param to compute power balance
flujoA(l, t)             param to compute power flow
fluMaxA(l, t)            param to check if max power flow is violated
fluMinA(l, t)            param to check if min power flow is violated
maxGenA(t, phiHT)        param to check if max power generation is violated
minGenA(t, phiHT)        param to check if min power generation is violated
chkBal(t, phiHT)         param to check if power balance is violated
RampDoA(t, phiHT)        param to check if ramp up is violated
RampUpA(t, phiHT)        param to check if ramp down is violated
infactible               aggregation of the num. of violated constraint /0/
countInfac(sim)          param. to storage the num. of violated cons. per sim
poe(sim)                 max percentage of deviation w.r.t. central value
flujoCont(k,l,t)         param. to compute power flow under contigency
flujoContMax(k,l,t)      param to check if max power flow under N-1 is violated
flujoContMin(k,l,t)      param to check if max power flow under N-1 is violated
;

table eWind(t, phiW)     deviation for wind generators
table eSolarg(t, phiS)   deviation for PV generators
table eLoad(t, i)        deviation for load
;

countInfac(sim) = eps;

loop (sim,
         bound = ord(sim)/card(sim)*2;
         poe(sim) = bound;

         loop(scen,
** create a random scenario
                 eWind(t, phiW) = bound*uniform(-1, 1)*wOut(t, phiW)*0.1;
                 eSolarg(t, phiS) = bound*uniform(-1, 1)*sOut(t, phiS)*0.1;
                 eLoad(t, i) = bound*uniform(-1, 1)*d(t, i)*0.05;

** compute the new power generation value according to the materializated scenario
                 pgA(t, phiHT) =
                         [p(t, phiHT)
                         + sum((phiW, tt)${ord(tt) le ord(t)}, gammaWind(phiHT, t, phiW, tt)*eWind(tt, phiW))
                         + sum((phiS, tt)${ord(tt) le ord(t)}, gammaSolar(phiHT, t, phiS, tt)*eSolarg(tt, phiS))
                         + sum((ii, tt)${ord(tt) le ord(t)}, gammaLoad(phiHT, t, tt)*eLoad(tt, ii))]*x(t, phiHT)
                 ;
** compute power balance
                 powBal(t) =
                         + sum(phiHHTT, pgA(t, phiHHTT))
                         + sum(phiW, wOut(t, phiW) + eWind(t, phiW))
                         + sum(phiS, sOut(t, phiS) + eSolarg(t, phiS))
                         - sum(i, d(t, i) + eLoad(t, i) - R(t, i))
                 ;
** compute power flow
                 flujoA(l, t) =
                         sum(i,
                                 PTDF(l, i)*[
                                 sum(phiHT, gen_map(phiHT, i) * pgA(t, phiHT))
                                 + sum[phiW, w_map(phiW, i) * (wOut(t, phiW) + eWind(t, phiW))]
                                 + sum[phiS, s_map(phiS, i) * (sOut(t, phiS) + eSolarg(t, phiS))]
                                 - d(t, i) - eLoad(t, i)
                                 + R(t, i)
                                 ]
                         )
               ;
** compute power flow under contingecies
                 flujoCont(k, l, t) =
                         sum(i,
                                OTDF(k, l, i)*[
                                  sum(phiHT, gen_map(phiHT, i) * pgA(t, phiHT))
                                + sum[phiW, w_map(phiW, i) * (wOut(t, phiW) + eWind(t, phiW))]
                                + sum[phiS, s_map(phiS, i) * (sOut(t, phiS) + eSolarg(t, phiS))]
                                - d(t, i) - eLoad(t, i)
                                + R(t, i)
                                ]
                        )
                 ;
** check violated constraints
                 flujoContMax(k, l, t) = 1${  flujoCont(k, l, t) gt  lineData(l, 'limit')*lineCap + tol};
                 flujoContMin(k, l, t) = 1${  flujoCont(k, l, t) lt -lineData(l, 'limit')*lineCap - tol};

                 fluMaxA(l, t) = 1${  flujoA(l, t) gt  lineData(l, 'limit')*lineCap + tol};
                 fluMinA(l, t) = 1${  flujoA(l, t) lt -lineData(l, 'limit')*lineCap - tol};

                 maxGenA(t, phiHT)${genData(phiHT, 'bus') ne idSlack} = 1${pgA(t, phiHT) gt genData(phiHT, 'Pmax')*x(t, phiHT) + tol};
                 minGenA(t, phiHT)${genData(phiHT, 'bus') ne idSlack} = 1${pgA(t, phiHT) lt genData(phiHT, 'Pmin')*x(t, phiHT) - tol};

                 chkBal(t, phiHT)${genData(phiHT, 'bus') eq idSlack} = 1${ powBal(t) lt  0 - tol };

                 RampDoA(t, phiHT) = 1${[ -genData(phiHT, 'Rdo') gt pgA(t, phiHT) - pgA(t-1, phiHT) + tol] and [ord(t) gt 1]};
                 RampUpA(t, phiHT) = 1${[pgA(t, phiHT) - pgA(t-1, phiHT) gt genData(phiHT, 'Rup') + tol] and [ord(t) gt 1]};
** agreggate all the violated constraint
                 infactible =
                 sum([t, phiHT],
                         + RampDoA(t, phiHT)
                         + RampUpA(t, phiHT)
                         + chkBal(t, phiHT)
                         + maxGenA(t, phiHT)
                         + minGenA(t, phiHT)
                 )
                 + sum([l, t],
                         fluMaxA(l, t) + fluMinA(l, t)
                         +sum(k, flujoContMax(k, l, t) + flujoContMin(k, l, t))
                 );

                 countInfac(sim) = countInfac(sim) + infactible;
         );
);

execute_unload 'MC_Past%past%.gdx';
